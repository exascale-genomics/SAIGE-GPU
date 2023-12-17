#!/bin/bash

#SBATCH -A med112
#SBATCH -N 50
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/logs/finemapping.downsampled_fine-map.EUR.parallel.run1_sorted.Anurag.out
#SBATCH -e /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/logs/finemapping.downsampled_fine-map.EUR.parallel.run1_sorted.Anurag.err
#SBATCH --distribution=cyclic:cyclic

module load python/3.7-anaconda3
module load r/4.0.3-py3

#################################################### Set global variables and file locations ##############################################################

#Set global variables
down_anc=EUR
ref_anc=AFR 

#Script locations
wrapper_script=/ccs/home/mconery/Pan-Trait_Colocalization/3_finemap/check_and_map_loci.py 
execution_script=/ccs/home/mconery/Pan-Trait_Colocalization/2_matrix_making/execute_command.py

#Files and directories
out_dir=/gpfs/alpine/proj-shared/med112/results/DOWNSAMPLING
locus_file=$out_dir/downsampling_loci.grch37.high_thresh.broad.csv
matrix_dir=$out_dir/matrix_dir_"$down_anc"_to_"$ref_anc" 
susie_dir=$out_dir/susie_dir_"$down_anc"_to_"$ref_anc" #Output directory
stats_dir=/gpfs/alpine/proj-shared/med112/task0101113/output/GIA_ANC_Run/downsample_summary_stats/QCed_Files #Directory of all summary stats
command_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/fine-map.downsample."$down_anc"_to_"$ref_anc".cmds.txt 
stats_name_conv="TRAIT.ANC.downsample.summary.txt"

###########################################################################################################################################################

#Make susie directory if it doesn't exist
mkdir -p $susie_dir
mkdir -p $susie_dir".2"

#Call python wrapper script with given parameters
python3 $wrapper_script -l $locus_file -c $command_file -o $susie_dir -m $matrix_dir -s $stats_dir -a $down_anc -n $stats_name_conv

#Iterate over the command file and run the lines
declare -a cmd_array=()
while read line; do
		cmd_array+=("$line")
done < $command_file

echo "JOB SIZE:"
echo "${#cmd_array[@]}"
for (( i=0; i<${#cmd_array[@]}; i++ )); do
	filename=${cmd_array[$i]}
        filename=${filename%___-o___*}
	filename=${filename#*___-p___}
	echo jobname=$i for $filename
	srun -n1 -c32 --exclusive --job-name=$i python3 $execution_script ${cmd_array[$i]} &
	if [ $i -ge 480 ]; then # (32/4)*60
	    wait
        fi
done

wait
