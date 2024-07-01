This is a Python Implementation of DirecAp Algorithm

To run this code use Python version >3.10 and the following command
python coexpression_networks_reconstruction.py ARG1 ARG2 ARG3 ARG4 ARG5 ARG6 ARG6 ARG7 ARG8 ARG9
with
ARG1: path to the original aracne jar file (provided in this repository).
ARG2: path to the modified aracne jar file (provided in this repository).
ARG3: Input file in the form of a tab delimited txt file with one row for each gene/protein and one column for each sample. Missing values are not allowed, while the input files's first column and row include the gene/protein ids and sample ids, respectively. 
ARG4: name of the folder that the results will be stored on.
ARG5: Tab delimited file with one column with protein ids that can be used as source in the reconstructed regulatory network.
ARG6: path to the folder where SIREN code has been placed (provided in this repository)
ARG7: Percentage of the average correlations of each molecule to use for the dynamic threshold calculation, default value: 0.1
ARG8: p-value threshold, default: 0.05
ARG9: Number of bootstraps to be used, default 100
