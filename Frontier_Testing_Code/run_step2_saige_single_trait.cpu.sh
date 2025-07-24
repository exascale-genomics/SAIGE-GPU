#!/bin/bash
# Parse arguments
bgen_file=$1
var_list=$2
sample_file=$3
model_file=$4
variance_ratio_file=$5
chrom=$6
output_prefix=$7
out_dir=$8

# Start timing
start_time=$(date +%s.%N)

#Copy over the files to the burst buffer
inp_dir=$(dirname $sample_file)
cp $bgen_file $inp_dir
cp $bgen_file.bgi $inp_dir
cp $var_list $inp_dir
cp $model_file $inp_dir
cp $variance_ratio_file $inp_dir

#Change over into output directory
cd $out_dir

# Path to Singularity image and SAIGE Step 2 script
saige_sif="/lustre/orion/bif154/proj-shared/mconery/tools/saige_1.4.4.1.sif"
saige_script="step2_SPAtests.R"

# Output file for this chromosome
saige_output="${output_prefix}.chr${chrom}.step2.txt"

# Optional: Set additional parameters
minMAF=0
minMAC=0.5
LOCO=TRUE
is_Firth_beta=TRUE
allele_order="ref-first"
is_fastTest=TRUE
pCutoffforFirth=0.05
n_Threads=$(($SLURM_CPUS_ON_NODE  / 2))

echo "NOTE: RUNNING on $n_Threads cores"

# Run SAIGE Step 2 with Singularity
singularity exec --bind $inp_dir:/input $saige_sif $saige_script \
    --bgenFile=/input/$(basename $bgen_file) \
    --bgenFileIndex=/input/$(basename $bgen_file).bgi \
    --sampleFile=/input/$(basename $sample_file) \
    --GMMATmodelFile=/input/$(basename $model_file) \
    --varianceRatioFile=/input/$(basename $variance_ratio_file) \
    --SAIGEOutputFile=$saige_output \
    --chrom=$chrom \
    --AlleleOrder=$allele_order \
    --minMAF=$minMAF \
    --minMAC=$minMAC \
    --is_Firth_beta=$is_Firth_beta \
    --pCutoffforFirth=$pCutoffforFirth \
    --LOCO=$LOCO \
    --idstoIncludeFile=/input/$(basename $var_list) \
    --nThreads=$n_Threads

#Clean temp files
rm $output_prefix.chr"$chrom".step2.txt[0-9]*

echo "SUCCESS: SAIGE step 2 completed for $saige_output"

# End timing and calculate duration
end_time=$(date +%s.%N)
duration=$(echo "$end_time - $start_time" | bc -l)

# Save timing results
echo "${output_prefix%.*},${output_prefix#*.},chr"$chrom",$duration" >> $out_dir/timing2.csv