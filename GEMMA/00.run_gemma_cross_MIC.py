"""
Command line tool for executing GEMMA

Authors:
Anna G. Green
"""
import pandas as pd
import numpy as np
import ruamel.yaml as yaml
import sys
import os
from copy import deepcopy
from codebase import *

config_file = sys.argv[1]

# Read in the config file and prepare
kwargs = yaml.safe_load(open(config_file, "r"))
gemma_executable = kwargs["gemma"]
drug = kwargs["drug"]
drug_abbrev = kwargs["drug_abbrev"]

# With gemma, anything after the -o (output) flag with automatically be prepended with "output/"
# Therefore, for saving in GEMMA with the -o flag use output_prefix, and for referring to those files, use prefix
prefix = f"output/{drug}_potentiators"
output_prefix = f"{drug}_potentiators"
os.system(f"mkdir {prefix}")
os.system(f"mkdir {prefix}/0.input_loci_files")
os.system(f"mkdir {prefix}/1.job_submission_files")
os.system(f"mkdir {prefix}/2.run_output")

#### Read in the genotype matrix and slice for isolates with relevant phenotype ####
isolates = pd.read_pickle(kwargs['isolate_annotation_file'])
print("reading the genotype matrix")
geno_mat = np.load(kwargs['genotype_array_file'])

# get the isolates with phenotypes to keep
phenotypes = pd.read_csv(kwargs['phenotype_file'], index_col=0)
phenotypes_found_all = phenotypes.loc[~phenotypes[f"{drug_abbrev}_midpoint"].isnull(),:]
phenotypes_found = phenotypes_found_all[["ROLLINGDB_ID", f"{drug_abbrev}_midpoint"]]
phenotypes_found[drug] = np.log(phenotypes_found[f"{drug_abbrev}_midpoint"])
isolates_with_phenotype = isolates.query("isolate_ID in @phenotypes_found.ROLLINGDB_ID")

print(f"found {len(isolates_with_phenotype)} isolates with phenotype for {drug}")

#### Prepare the phenotype input file ####
print("writing the phenotype file")
pheno_file = f"{prefix}/{drug}.phenotypes"
write_phenotype_file_MIC(pheno_file, phenotypes_found, isolates_with_phenotype, drug)

#### If the .CXX file does not exist, create it ####
cxx_file = f"{prefix}/MAF_0p001.cXX.txt"
U_file = f"{prefix}/MAF_0p001.eigenU.txt"
relatedness_loci_file = f"{prefix}/0.input_loci_files/MAF_0p001.loci"

if not os.path.isfile(cxx_file) or not os.path.isfile(U_file):

    print("generating the CXX file", relatedness_loci_file)
    # get the SNPs with 0.1% MAF for constructing the locus file
    snps = pd.read_csv(kwargs['snp_annotation_file'], index_col=0)

    # slice the matrix for only isolates and snps of interest
    matrix = geno_mat[snps.index, :][:, isolates_with_phenotype.index]
    matrix = translate_genotypes(matrix, genotype_translator).T
    print("Number of gaps in the matrix", np.sum(matrix=='-'))


    # Add the ancestor sequence
    matrix, isolates = add_ancestor(kwargs["mtb_ancestor_fasta"], snps, matrix, isolates_with_phenotype)

    print(" the shape of the input matrix", matrix.shape)
    print("number of isolates", len(isolates_with_phenotype))
    # Generate the locus file
    prepare_loci_file(matrix, isolates, list(snps.pos), relatedness_loci_file, prefix=prefix, gap_code="NA")

    # Generate the cXX file
    outfile = f"{output_prefix}/MAF_0p001"
    string1 = f"{gemma_executable} -g {relatedness_loci_file} -p {pheno_file} -gk -maf 0 -o {outfile}"
    print(string1)

    # Run the eigen decomposition
    outfile =f'{output_prefix}/MAF_0p001'
    string2 = f"{gemma_executable} -g {relatedness_loci_file} -p {pheno_file} -k {cxx_file} -maf 0 -eigen -notsnp -o {outfile}"
    print(string2)

    print("creating GRM")
    os.system(string1)

    print("running eigendecomp")
    os.system(string2)

## Prepare locus pair cross runs ###

# Read in table of all pairs that need testing
df = pd.DataFrame()
for file in kwargs["input_pairs_file"]:
    df = pd.concat([df,pd.read_csv(file)])
    
print(f"Submitting cross for {len(df)} pairs")

complete_snps = pd.read_pickle("../input/genotypes_SNP_annotation.pkl")

table = [] # i, j, n_i, n_j, n_ij

# iterates through all pairs and executes a GEMMA run for that pair using sbatch
for _, row in df.iterrows():
    position_i = row.position_i
    position_j = row.position_j

    print("setting up cross for ", position_i, position_j)

    # Need to query individually and then concatenate to preserve i,j order
    snps_subset = pd.concat([complete_snps.query("pos==@position_i"),complete_snps.query("pos==@position_j")])

    matrix = geno_mat[snps_subset.index, :][:, isolates_with_phenotype.index]
    matrix = translate_genotypes(matrix, genotype_translator).T

    # add the ancestral sequence
    matrix, isolates = add_ancestor(kwargs["mtb_ancestor_fasta"], snps_subset, matrix, isolates_with_phenotype)

    # compute the major allele matrix
    major_allele_matrix = reference_allele_matrix_floats(matrix,np.array(matrix[0,:]))

    # compute the multiplicative term nd add to the matrix
    and_term = np.multiply(major_allele_matrix[:,0], major_allele_matrix[:,1]).reshape(-1,1)
    run_matrix = np.concatenate([major_allele_matrix, and_term], axis=1)

    table.append([position_i, position_j, np.nansum(major_allele_matrix[:,0]), np.nansum(major_allele_matrix[:,1]), np.nansum(and_term)])
    # prepare and save the locus file
    positions = [position_i, position_j, f"{position_i}_{position_j}"]
    trait_loci_file =f'{prefix}/0.input_loci_files/{position_i}_{position_j}_cross.loci'
    prepare_loci_file(run_matrix, isolates, positions, trait_loci_file, matrix_mode="two_coded")

    # run the actual gwas
    outfile =f'{output_prefix}/2.run_output/{position_i}_{position_j}_cross'
    string = f"{gemma_executable} -g {trait_loci_file} -p {pheno_file} -d {prefix}/MAF_0p001.eigenD.txt -u {prefix}/MAF_0p001.eigenU.txt -maf 0 -miss 0.2 -lmm -notsnp -o {outfile}"
    job_file = f'{prefix}/1.job_submission_files/{position_i}_{position_j}_cross.sh'
    stdout_str = f'{prefix}/2.run_output/{position_i}_{position_j}_cross'
    string_to_write = job_submission_string.format(prefix=stdout_str, command=string)

    with open(job_file, "w") as out:
        out.write(string_to_write)

    os.system(f"sbatch {job_file}")


pd.DataFrame(table, columns=["position_i", "position_j", "N_i", "N_j", "N_ij"]).to_csv(f"{prefix}/{drug}_mutation_count_table.csv")
