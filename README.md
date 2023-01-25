# DependentMutations

This repository contains the series of Jupyter notebooks used to generate the results for Green et al, 2021, forthcoming. 

The notebooks are run in numerical order. 
  - **00.convert_tree_data_to_event_matrices.ipynb** - converts the SNPPar output data into matrices of L x N, where L is the number of branches and N is the number of SNPs. 
  Computes a matrix of mutation events, reversions events, and genetic background, for each lineage individually and then combines. 

  - **01.event_matrices_to_pairs.ipynb** - Computes the first-order (single site) number of mutations. Computes matrix of all possible pairs of mutations. 
    Determines which pairs of mutations are found to ever co-occur. 
  
  - **02.compute_p_values.ipynb** - For pairs of mutations that actually co-occurr, computes the p-value and Benjamini-Hochberg multiple testing correction. 
  
  - **03.analyze_pvalues.ipynb** - Creates tables of hits that are related to antibiotic resistance


### Download the intermediate and output data here: <future_link>

### Download the input data here: <future_link>

Note that the input data is from two other publications that should be cited accordingly: 

  - Comas et al, Nat Gen, 2013 10.1038/ng.1038 (M. tuberculosis pan-susceptible ancestor sequence)

  - Vargas et al, 2021, forthcoming (phylogenies and SNPpar output ancestral sequence reconstructions)
