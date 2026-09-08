These scripts were created to help speed up and automate next-generation sequencing data processing.

Each step of this pipeline will skip if the 1st sample's 1st replicate's associated output file exists. If the previous run of the script failed in that step, please delete at least
the 1st sample's 1st replicate associated output file. 

These scripts script expects the user to have a conda environment available named "deeptools_env" with the deepTools suite installed in said environment.
