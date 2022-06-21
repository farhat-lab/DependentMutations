"""
Command line tool for executing GEMMA

Authors:
Anna G. Green
"""
import pandas as pd
import numpy as np
import ruamel.yaml as yaml
import sys
from evcouplings.utils import valid_file, verify_resources, run
from evcouplings.align import Alignment
import os
from copy import deepcopy
from scipy.stats import mode

ALPHABET_NOGAP = "ACGT"
GAP_CHAR = "-"
ALPHABET = GAP_CHAR + ALPHABET_NOGAP

genotype_translator = {
    0:"A", 1:"C", 2:"G", 3:"T", 9:"-"
}

job_submission_string = "#!/bin/bash\n"+\
"#SBATCH -c 1                               # Request one core\n"+\
"#SBATCH -t 0-11:05                         # Runtime in D-HH:MM format\n"+\
"#SBATCH -p short                           # Partition to run in\n"+\
"#SBATCH --mem=80G                         # Memory total in MB (for all cores)\n"+\
"#SBATCH -o {prefix}_%j.out                 # File to which STDOUT will be written, including job ID (%j)\n"+\
"#SBATCH -e {prefix}_%j.err                 # File to which STDERR will be written, including job ID (%j)\n"+\
"{command}"


def map_matrix(matrix, map_):
    """
    Map elements in a numpy array using alphabet
    Parameters
    ----------
    matrix : np.array
        Matrix that should be remapped
    map_ : defaultdict
        Map that will be applied to matrix elements
    Returns
    -------
    np.array
        Remapped matrix
    """
    return np.vectorize(map_.__getitem__)(matrix)


def calculate_frequencies(matrix, num_symbols):
    """
    Calculate single-site frequencies of symbols in alignment
    Parameters
    ----------
    matrix : np.array
        N x L matrix containing N sequences of length L.
        Matrix must be mapped to range(0, num_symbols) using
        map_matrix function
    num_symbols : int
        Number of different symbols contained in alignment
    Returns
    -------
    np.array
        Matrix of size L x num_symbols containing relative
        column frequencies of all characters
    """
    N, L = matrix.shape
    fi = np.zeros((L, num_symbols))
    for s in range(N):
        for i in range(L):
            fi[i, matrix[s, i]] += 1
    return np.divide(fi, N)

def write_phenotype_file_MIC(gemma_input_phenotype_file, phenotypes_found, isolates_with_phenotype, drug):

    """
    Convert input phenotype file into correct format for GEMMA and saves file.
    Adds an ancestral sequence and assumes it is sensitive

    Parameters
    ---------
    gemma_input_phenotype_file: str
        name of file to save

    phenotypes_found: pd.DataFrame

    isolates: pd.DataFrame

    """

    phenotypes_found_subset = phenotypes_found.drop_duplicates(subset=["ROLLINGDB_ID"])
    phenotypes = isolates_with_phenotype[["isolate_ID"]].merge(
        phenotypes_found_subset[["ROLLINGDB_ID", drug]], right_on="ROLLINGDB_ID", left_on="isolate_ID", how="left"
    )

    # Add the ancestor sequence, assuming ancestor is sensitive
    ancestor = pd.DataFrame({
        "Isolate_ID": ["MT_H37Rv_anc"],
        "ROLLINGDB_ID": [None],
        drug: phenotypes[drug].min()
    })
    phenotypes = pd.concat([ancestor, phenotypes])

    phenotypes.loc[:,drug].to_csv(gemma_input_phenotype_file, header=False, index=None, na_rep="NA")

def calculate_reference_allele_matrix(matrix, reference_alleles, encode_gaps=True, gap_code=np.nan,
        ref_allele_code=0, non_ref_allele_code=1):
    """
    Encodes a mapped genotype matrix with reference and non-reference allele codes
    NOTE: This function currently only works for 5-character (with gaps) alphabet

    ----------
    matrix : np.array
        N x L matrix containing N sequences of length L.
        Matrix must be mapped to range(0, num_symbols) using
        map_matrix function
    reference_alleles : np.array
        1 x L matrix with coding of reference allele for each position in L
    encode_gaps: boolean, optional (default=True)
        If True, gaps will be fille with gap_code. If False, gaps will be
        eliminated (ie, replaced with the major allele code)
    gap_code: str or numeric, optional (default=np.nan)
        encoding to use for gap characters if nogap=False
    ref_allele_code: str or numeric, optional (default=0)
        encoding to use for major alleles
    non_ref_allele_code: str or numeric, optional (default=1)
        encoding to use for minor alleles

    Returns
    -------
    np.array
        Matrix of size N x L with reference/non-reference allele encoding
        for all positions in L
    """

    if encode_gaps:
        gap_code=gap_code
    else:
        gap_code=ref_allele_code

    encoded_matrix = np.zeros(shape=matrix.shape, dtype=object)

    # For each column
    for i in np.arange(0, encoded_matrix.shape[1]):

        #grab that column
        vec = matrix[:, i]

        # create an encoded version of the column, with minor allele code as default
        encoded_vec = np.array([non_ref_allele_code] * len(vec), dtype=object)

        # If the row has a major allele, replace in encoded vec
        is_ref_allele = np.equal(vec, reference_alleles[i])
        encoded_vec[is_ref_allele] = ref_allele_code

        # If the row has a gap character, replace in encoded vec
        is_gap = np.isnan(vec.astype(float))
        encoded_vec[is_gap] = gap_code

        # save in the encoded matrix
        encoded_matrix[:, i] = encoded_vec

    return encoded_matrix

class GenotypeMatrix:
    """
    Takes inspiration from the evcouplings python package Alignment class
    """
    def __init__(self, strain_list, position_list, genotype_data=None, matrix_data=None):

        self.positions = position_list
        self.position_to_index = {p:idx for idx, p in enumerate(self.positions)}
        self.strains = strain_list
        self.strain_to_index = {p:idx for idx, p in enumerate(self.strains)}

        # If user provided genotype data, do some sanity checks
        self.genotype_data = genotype_data
        # Make sure that every position and every genotype is represented in the genotype_data
        if not self.genotype_data is None:
            assert set(self.positions).issubset(genotype_data.POS)
            assert set(self.strains).issubset(genotype_data.STRAIN)

        # if user provided matrix data, set it up
        if not matrix_data is None:
            assert matrix_data.shape[0] == len(self.strains)
            assert matrix_data.shape[1] == len(self.positions)
            self.matrix = matrix_data

            self.reference_alleles = deepcopy(self.matrix[0, :])

            if 0 in matrix_data:
                self.matrix_mapped = matrix_data
            else:
                self.matrix_mapped = None

        else:
            self.matrix = np.zeros((len(self.strains), len(self.positions)), dtype=str)
            self.reference_alleles = None
            self.matrix_mapped = None

        self.alphabet_map = {
            c: i for i, c in enumerate(ALPHABET)
        }
        self.num_symbols = len(self.alphabet_map)

        self._frequencies = None
        self._major_alleles = None


    @classmethod
    def from_df(cls, file):
        """
        Initializes a GenotypeMatrix using an input pd.DataFrame, where rows = strains and columns=positions, and
        the data will become the contents of self.matrix
        """
        df = pd.read_csv(file, index_col=0, header=0)

        return cls(
            list(df.index), list(df.columns), genotype_data=None, matrix_data=df.values
        )

    def make_matrix(self):
        '''
        Populates the matrix of strains by positions with genotype data from
        the input genotype data file via the following heuristic:

        0. Initialize every position in every genome as the reference allele
        1. For each position in each strain, determine if there is any information
            in the genotype data input. If not, leave that position as Reference. If yes,
            proceed to step 2.
        2. If the alternate allele is an insertion or a non-standard character, insert a gap
            at that position.
        3. If the alternate allele did not pass quality threshold (determined by QUAL flag being 0),
            insert a gap at that position
        4. If alternate allele is a SNP and is not low quality, replace reference with alternate allale

        Returns
        -------
        np.ndarray:
            an N_strains x N_positions array of strings
        '''
        #print("preparing the reference")

        position_subset = self.genotype_data.query("POS in @self.positions")

        position_groups = position_subset.groupby("POS")
        for position, subset in position_groups:
            allele = subset.loc[subset.index[0], "REF"]
            self.matrix[:, self.position_to_index[position]] = allele

        self.reference_alleles = deepcopy(self.matrix[0, :])
        #print("assigned the reference", self.reference_alleles)

        # iterate through each strain
        strain_group = position_subset.query("STRAIN in @self.strains").groupby("STRAIN")

        for strain, subset in strain_group:

            # iterate through each position
            position_groups = subset.groupby("POS")
            for position, subsubset in position_groups:

                # if we didn't find any information for this position, leave as reference allele
                if len(subsubset) == 0:
                    continue

                subsubset = subsubset.iloc[0,:]

                # If we did find information but its an indel, treat as missing
                if len(subsubset.ALT) > 1 or  not subsubset.ALT in ALPHABET:
                    self.matrix[self.strain_to_index[strain],self.position_to_index[position]] = "-"

                # If we found information but it didn't pass quality filter, treat as missing
                elif subsubset.QUAL==0:
                    self.matrix[self.strain_to_index[strain],self.position_to_index[position]] = "-"

                # else, change to alterante allele
                else:
                    self.matrix[self.strain_to_index[strain],self.position_to_index[position]] = subsubset.ALT


    def __ensure_mapped_matrix(self):
        """
        Ensure self.matrix_mapped exists
        """
        if self.matrix_mapped is None:
            self.matrix_mapped = map_matrix(
                self.matrix, self.alphabet_map
            )

    @property
    def major_alleles(self):
        """
        Returns/calculates identity of major allele for each position in alignment.
        Also sets self._major_alleles member variable for later reuse.

        Returns
        -------
        np.array
            Reference to self._major_alleles
        """
        if self._major_alleles is None:

            self._major_alleles = np.argmax(
                self.frequencies, axis=1
            )

        return self._major_alleles

    @property
    def frequencies(self):
        """
        Returns/calculates single-site frequencies of symbols in alignment.
        Also sets self._frequencies member variable for later reuse.

        Returns
        -------
        np.array
            Reference to self._frequencies
        """
        if self._frequencies is None:
            self.__ensure_mapped_matrix()

            self._frequencies = calculate_frequencies(
                self.matrix_mapped, self.num_symbols
            )

        return self._frequencies

    def write_fasta(self, fileobj, strain_list=None):
        """
        Creates a fasta-formatted file of pseudo-sequences for strains and positions in self.matrix

        Parameters
        ----------
        fileobj: filehandle
            fileobject to which to write
        strain_list: list of str, optional (default=None)
            List of strain names, in order, which will be written to file. If None,
            will write all strains in self.strains

        """
        if strain_list is None:
            strain_list = self.strains

        for strain in strain_list:
            seq = "".join(self.matrix[self.strain_to_index[strain],:])
            fileobj.write(f">{strain}\n{seq}\n")

    def write_GEMMA(self, fileobj, positions=None, matrix_mode="four_encoded", major_allele_string="\"A\"", minor_allele_string="TRUE",
        major_allele_code='0', minor_allele_code='1',gap_code="NA"):
        """
        Creates a file of loci in correct input format for GEMMA

        User warning: Changing the optional variables to values other than their defaults may result in errors when
        using the program GEMMA. Verify that all inputs are still valid before proceeding.

        Gemma input format:
        position_name,major_allele_string,minor_allele_string,{allele code for strain 1, ..., allele code for strain N}

        Parameters
        ----------
        fileobj: filehandle
            file to which to write
        positions: list of str, optional (default: None)
            specify subset of positions to write. If None, use all positions
        major_allele_string: str, optional (default "A")
            string used to represent the major allele
        minor_allele_string: str, optional (default "T")
            string used to represent the minor allele
        major_allele_code: str, optional (default "1")
            string used to code for the major allele in the input
        minor_allele_code: str, optional (default "0")
            string used to code for the major allele in the input
        gap_code: str, optional (default "NaN")
            string used to code for the gap/missing data in the input

        """

        if positions is None:
            positions = self.positions

        self.__ensure_mapped_matrix()
        if not 0 in self.reference_alleles:
            reference_alleles = [self.alphabet_map[x] for x in self.reference_alleles]
        else:
            reference_alleles = self.reference_alleles

        if matrix_mode=='four_encoded':
            major_allele_matrix = calculate_reference_allele_matrix_for_pair_regression(
                    self.matrix_mapped,
                    reference_alleles,
                    encode_gaps=True,
                    gap_code=gap_code,
                    ref_allele_code=major_allele_code,
                    non_ref_allele_code=minor_allele_code
                )
        else:
            print("using two-encoded matrix mode")
            major_allele_matrix = float_matrix_to_object(
                    self.matrix_mapped,
                    reference_alleles,
                    encode_gaps=True,
                    gap_code=gap_code,
                    ref_allele_code=major_allele_code,
                    non_ref_allele_code=minor_allele_code
                )

        for position in positions:
            # slice corresponding to all strains
            vec = major_allele_matrix[:, self.position_to_index[position]]

            fileobj.write(f'{position},{major_allele_string},{minor_allele_string},{",".join(vec)}\n')

    def to_df(self):
        """
        Creates a pd.DataFrame, with rows corresponding to strains, columns corresponding to positions,
        and entries corresponding to the alleles in self.matrix
        Returns
        -------
        pd.DataFrame
        """
        df = pd.DataFrame(self.matrix, index=self.strains, columns=self.positions)
        return df

def translate_genotypes(matrix, translator):
    """
    Translate matrix of genotypes into string characters. Used for converting Roger's numeric genotype
    matrix back into strings of ATCG-

    Parameters
    ----------

    matrix: np.array
        input genotype matrix
    translator: dict
        translator for converting entires of matrix to strings
    Returns
    -------
    np.array of str, translated matrix

    """

    translated_matrix = np.zeros_like(matrix, dtype='str')
    for i in range(matrix.shape[0]):
        translated_matrix[i, :] = [translator[x] for x in matrix[i,:]]
    return translated_matrix


# Add the ancestral sequence to the sequence matrix
def add_ancestor(mtb_ancestor_fasta, snps_to_keep, matrix, isolates_with_phenotype):
    mtb_ancestor = Alignment.from_file(open(mtb_ancestor_fasta, "r"))
    indices_to_slice = [x-1 for x in snps_to_keep.pos]
    ancestor_sequence = mtb_ancestor.matrix[:, indices_to_slice]

    matrix  = np.concatenate([ancestor_sequence, matrix])
    isolate_labels = ["MTB_ancestor"] + list(isolates_with_phenotype.isolate_ID)

    return matrix, isolate_labels

def write_phenotype_file(gemma_input_phenotype_file, phenotypes_found, isolates_with_phenotype, drug):

    """
    Convert input phenotype file into correct format for GEMMA and saves file.
    Adds an ancestral sequence and assumes it is sensitive

    Parameters
    ---------
    gemma_input_phenotype_file: str
        name of file to save

    phenotypes_found: pd.DataFrame

    isolates: pd.DataFrame

    """

    phenotypes = isolates_with_phenotype[["isolate_ID"]].merge(
        phenotypes_found[["Isolate", drug]], left_on="isolate_ID", right_on="Isolate", how="left"
    )

    # Add the ancestor sequence, assuming ancestor is sensitive
    ancestor = pd.DataFrame({
        "Isolate": ["MT_H37Rv_anc"],
        "isolate_ID": [None],
        drug: "S"
    })
    phenotypes = pd.concat([ancestor, phenotypes])

    encoding = {"R": 1, "S": 0}
    phenotypes[drug] = [encoding[x] for x in phenotypes[drug]]

    phenotypes.loc[:,drug].to_csv(gemma_input_phenotype_file, header=False, index=None, na_rep="NA")

def calculate_reference_allele_matrix_for_pair_regression(matrix, reference_alleles, encode_gaps=True, gap_code=np.nan,
        ref_allele_code=0, non_ref_allele_code=1):
    """
    Encodes a mapped genotype matrix with reference and non-reference allele codes
    NOTE: This function currently only works for 5-character (with gaps) alphabet

    ----------
    matrix : np.array
        N x L matrix containing N sequences of length L.
        Matrix must be mapped to range(0, num_symbols) using
        map_matrix function
    reference_alleles : np.array
        1 x L matrix with coding of reference allele for each position in L
    encode_gaps: boolean, optional (default=True)
        If True, gaps will be fille with gap_code. If False, gaps will be
        eliminated (ie, replaced with the major allele code)
    gap_code: str or numeric, optional (default=np.nan)
        encoding to use for gap characters if nogap=False
    ref_allele_code: str or numeric, optional (default=0)
        encoding to use for major alleles
    non_ref_allele_code: str or numeric, optional (default=1)
        encoding to use for minor alleles

    Returns
    -------
    np.array
        Matrix of size N x L with reference/non-reference allele encoding
        for all positions in L
    """

    if encode_gaps:
        gap_code=gap_code
    else:
        gap_code=ref_allele_code

    encoded_matrix = np.zeros(shape=matrix.shape, dtype=object)

    # For each column
    for i in np.arange(0, encoded_matrix.shape[1]):

        #grab that column
        vec = matrix[:, i]

        # create an encoded version of the column, with minor allele code as default
        encoded_vec = np.array([non_ref_allele_code] * len(vec), dtype=object)

        # If the row has a major allele, replace in encoded vec

        is_ref_allele = vec == reference_alleles[i]
        encoded_vec[is_ref_allele] = ref_allele_code

        # If the row has a gap character, replace in encoded vec
        is_gap = vec == 0
        encoded_vec[is_gap] = gap_code

        # save in the encoded matrix
        encoded_matrix[:, i] = encoded_vec

    return encoded_matrix

def float_matrix_to_object(matrix, reference_alleles, gap_code="NA", encode_gaps=True,
        ref_allele_code=0, non_ref_allele_code=1):
    """
    Converts a 3-character float matrix into an object for writing

    0: input reference allele codes
    1: input alternate allele codes
    np.nan: input gap code
    """

    encoded_matrix = np.zeros(shape=matrix.shape, dtype=object)

    # For each column
    for i in np.arange(0, encoded_matrix.shape[1]):

        #grab that column
        vec = matrix[:, i]

        # create an encoded version of the column, with minor allele code as default
        encoded_vec = np.array([non_ref_allele_code] * len(vec), dtype=object)

        # If the row has a major allele, replace in encoded vec

        is_ref_allele = vec == reference_alleles[i]
        encoded_vec[is_ref_allele] = ref_allele_code

        # If the row has a gap character, replace in encoded vec
        is_gap = np.isnan(vec)
        encoded_vec[is_gap] = gap_code

        # save in the encoded matrix
        encoded_matrix[:, i] = encoded_vec

    return encoded_matrix

def reference_allele_matrix_floats(matrix, reference_alleles, encode_gaps=True, gap_code=np.nan,
        ref_allele_code=0, non_ref_allele_code=1):
    """
    Encodes a mapped genotype matrix with reference and non-reference allele codes
    NOTE: This function currently only works for 5-character (with gaps) alphabet

    ----------
    matrix : np.array
        N x L matrix containing N sequences of length L.
        Matrix must be mapped to range(0, num_symbols) using
        map_matrix function
    reference_alleles : np.array
        1 x L matrix with coding of reference allele for each position in L
    encode_gaps: boolean, optional (default=True)
        If True, gaps will be fille with gap_code. If False, gaps will be
        eliminated (ie, replaced with the major allele code)
    gap_code: str or numeric, optional (default=np.nan)
        encoding to use for gap characters if nogap=False
    ref_allele_code: str or numeric, optional (default=0)
        encoding to use for major alleles
    non_ref_allele_code: str or numeric, optional (default=1)
        encoding to use for minor alleles

    Returns
    -------
    np.array
        Matrix of size N x L with reference/non-reference allele encoding
        for all positions in L
    """

    if encode_gaps:
        gap_code=gap_code
    else:
        gap_code=ref_allele_code

    encoded_matrix = np.zeros(shape=matrix.shape, dtype=float)

    # For each column
    for i in np.arange(0, encoded_matrix.shape[1]):

        #grab that column
        vec = matrix[:, i]

        # create an encoded version of the column, with minor allele code as default
        encoded_vec = np.array([non_ref_allele_code] * len(vec), dtype=float)

        # If the row has a major allele, replace in encoded vec

        is_ref_allele = vec == reference_alleles[i]
        encoded_vec[is_ref_allele] = ref_allele_code

        # If the row has a gap character, replace in encoded vec
        is_gap = vec =='-'
        encoded_vec[is_gap] = gap_code

        # save in the encoded matrix
        encoded_matrix[:, i] = encoded_vec

    return encoded_matrix


def prepare_loci_file(matrix,
    isolates, snps,  output_loci_file, matrix_mode="four_encoded", gap_code="NA", prefix=None):
    """
    Perpares input file for GEMMA
    Saves major and minor allele files if prefix is provided

    Parameters
    ----------
    matrix: np.array
        array of genotypes, N isolates x N positions
    isolates: iterable
        list of isolate names
    snps: iterable
        list of position names
    output_loci_file: str, optional (default=None)
        path to save the output file to
    gap_code: str, optional (default=NA)
        code with which to encode gaps in the GEMMA input loci file
    prefix: str, optional (default=None)
        output prefix for major allele and strains file to save
    """

    gm = GenotypeMatrix(
        isolates,
        snps,
        matrix_data = matrix
    )

    # write the input loci file
    of = open(output_loci_file, "w")
    gm.write_GEMMA(of, matrix_mode=matrix_mode, gap_code=gap_code)


    if prefix:
        # write the major allele file
        major_allele_file = f"{prefix}_major_alleles.csv"
        maj_alleles = pd.DataFrame(gm.major_alleles, index=gm.positions, columns=["major_alleles"])
        maj_alleles.loc[:,"major_alleles"] = [ALPHABET[x] for x in maj_alleles.major_alleles]
        maj_alleles.to_csv(major_allele_file)

        # write the strain file
        strain_file= f"{prefix}_strains.csv"
        strains = pd.DataFrame(gm.strains, columns=["strain"])
        strains.to_csv(strain_file, index=None)
