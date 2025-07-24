#!/bin/bash
#Extract command line parameters
plink_file=$1
pheno_file=$2
phenotype=$3
out_dir=$4
trait_type=$5

# Start timing
start_time=$(date +%s.%N)

#Change into the desired output directory
cpu_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/CPU_test/AFR"
cd $cpu_dir
#Copy in the necessary genotype files
inp_dir=$(dirname $plink_file)
cp $pheno_file $inp_dir
#Change into the output directory
cd $out_dir
#Get basename of plink file
plink_name=$(basename $plink_file)

#Check for trait type and then execute step 1
if [ $trait_type == "binary" ]; then
	singularity exec --bind $inp_dir:/data /lustre/orion/bif154/proj-shared/mconery/tools/saige_1.4.4.1.sif step1_fitNULLGLMM.R --plinkFile=/data/$plink_name --phenoFile=/data/$phenotype.$trait_type.txt --invNormalize=FALSE --phenoCol=PHENO --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --sampleIDColinphenoFile=IID --traitType=$trait_type --outputPrefix=$phenotype.$trait_type.step1_output --minMAFforGRM=0.01 --LOCO=T --IsOverwriteVarianceRatioFile=TRUE
else
	singularity exec --bind $inp_dir:/data /lustre/orion/bif154/proj-shared/mconery/tools/saige_1.4.4.1.sif step1_fitNULLGLMM.R --plinkFile=/data/$plink_name --phenoFile=/data/$phenotype.$trait_type.txt --invNormalize=TRUE --phenoCol=PHENO --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --sampleIDColinphenoFile=IID --traitType=quantitative --outputPrefix=$phenotype.$trait_type.step1_output --minMAFforGRM=0.01 --LOCO=T --IsOverwriteVarianceRatioFile=TRUE
fi


# End timing and calculate duration
end_time=$(date +%s.%N)
duration=$(echo "$end_time - $start_time" | bc -l)

# Save timing results
echo "$phenotype,$trait_type,$duration" >> $out_dir/timing.csv