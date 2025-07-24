#!/bin/bash
#Extract command line parameters
plink_file=$1
pheno_file=$2
prefix=$3
out_dir=$4
path_to_saige=$5
trait_type=$6

# Start timing
start_time=$(date +%s.%N)

#Run step 1
if [ $trait_type == "binary" ]; then
	Rscript $path_to_saige/extdata/step1_fitNULLGLMM.R --plinkFile=$plink_file --phenoFile=$pheno_file --invNormalize=FALSE --phenoCol=PHENO --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --sampleIDColinphenoFile=IID --traitType=$trait_type --outputPrefix=$out_dir/$prefix.step1_output --minMAFforGRM 0.01 --LOCO T  --IsOverwriteVarianceRatioFile=TRUE --nThreads=1
else
	Rscript $path_to_saige/extdata/step1_fitNULLGLMM.R --plinkFile=$plink_file --phenoFile=$pheno_file --invNormalize=TRUE --phenoCol=PHENO --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --sampleIDColinphenoFile=IID --traitType=quantitative --outputPrefix=$out_dir/$prefix.step1_output --minMAFforGRM 0.01 --LOCO T  --IsOverwriteVarianceRatioFile=TRUE --nThreads=1
fi

# End timing and calculate duration
end_time=$(date +%s.%N)
duration=$(echo "$end_time - $start_time" | bc -l)

# Save timing results
echo "$prefix,$trait_type,$duration" >> $out_dir/timing.csv