# Fine-Mapping Analysis Scripts

This directory contains all the scripts used to complete the fine-mapping analysis in Verma et al. The scripts were run on the Oak Ridge Leadership Computing Facility's Andes supercomputer. The scripts are broken down into the following five steps:

## 1) Locus Defining ##
These scripts define the loci used for fine-mapping by tiling the genome and merging together adjacent tiles with significant variants. The key shell script in this directory, **define_sig_loci.sh**, will submit the locus defining-process and merge the loci together into a comma-delimited list upon completion.
