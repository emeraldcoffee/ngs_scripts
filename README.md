These scripts were created to help speed up and automate next-generation sequencing data processing.

Expected file structures are shown as a comment at the top of each script.

Each step of this pipeline will skip if the 1st sample's 1st replicate's associated output file exists. If the previous run of the script failed in that step, please delete at least
the 1st sample's 1st replicate associated output file. 

The DNA sequencing scripts expect the user to have a conda environment available named "deeptools_env" with the deepTools suite installed in said environment.

Note that if the user is using spike-in normalization for DNA sequencing, either spike.fa or spike_index (from a previous spike-in normalization run) should be available in the usr directory.
