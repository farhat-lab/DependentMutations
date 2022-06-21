# DependentMutations

This repository contains the series of Jupyter notebooks used to generate the results for Green et al, 2021, forthcoming. 

The notebooks are run in numerical order. 
  - **00.convert_tree_data_to_event_matrices.ipynb** - converts the SNPPar output data into matrices of L x N, where L is the number of branches and N is the number of SNPs. 
  Computes a matrix of mutation events, reversions events, and genetic background, for each lineage individually and then combines. 

  - **01.event_matrices_to_pairs.ipynb** - Computes the first-order (single site) number of mutations. Computes matrix of all possible pairs of mutations. 
    Determines which pairs of mutations are found to ever co-occur. 
  
  - **02A.compute_p_values_sequential.ipynb** - For pairs of mutations that actually co-occurr on sequential branches, computes the p-value and Benjamini-Hochberg multiple testing correction. 
  
  - **02A.compute_p_values_sequential.ipynb** - For pairs of mutations that actually co-occurr on simultaneous branches, computes the p-value and Benjamini-Hochberg multiple testing correction. 
  
  - **03.single_mutation_annotation.ipynb** - Annotates all homoplastic SNPs with information about function, lineage, etc. Calculates enrichment statistics for paper
  
  - **04A.analyze_pvalues_sequential.ipynb** - analysis of significant, sequentially ocurring dependent mutation pairs

  - **04B.analyze_pvalues_sequential.ipynb** - analysis of significant, simultaneously ocurring dependent mutation pairs

  - **05.analyze_antibiotic_hits.ipynb** - analyze significant dependent mutations pairs related to antibiotic resistance

  - **06.antibiotic_potentiators.ipynb** - determination of which hits "potentiate" antibiotic resistance

  - **07.antibiotic_potentiators_table.ipynb** - creates tables of hits that potentiate antibiotic resistance for input to GEMMA

### Running GEMMA analyses

### Citations

Note that some of the input data is from other publications that should be cited accordingly: 

  - Comas et al, Nat Gen, 2013 10.1038/ng.1038 (M. tuberculosis pan-susceptible ancestor sequence)

  - Vargas et al, 2022, https://doi.org/10.1101/2022.06.10.495637 (phylogenies and SNPpar output ancestral sequence reconstructions)

  - Freschi et al, 2021, https://doi.org/10.1038/s41467-021-26248-1 (lineage SNP barcode)

  - Coll et al, 2014, https://doi.org/10.1038/ncomms5812 (lineage SNP barcode)

  - Kapopoulou A, Lew JM, Cole ST, Tuberculosis (Edinb), 2011. 10.1016/j.tube.2010.09.006 (Gene location data)

  - Minato et al, 2019, 10.1128/mSystems.00070-19, (Gene essentiality data)

  - WHO catalog of resistance variants, https://apps.who.int/iris/handle/10665/341906
