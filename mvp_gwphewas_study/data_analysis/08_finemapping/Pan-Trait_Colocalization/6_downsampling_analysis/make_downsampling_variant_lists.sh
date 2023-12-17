#!/bin/bash

#SBATCH -A med112
#SBATCH -N 60
#SBATCH -J make_downsampling_variant_lists
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/logs/finemapping.make_downsampling_variant_lists.parallel.out
#SBATCH -e /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/logs/finemapping.make_downsampling_variant_lists.parallel.err
#SBATCH --distribution=cyclic:cyclic

module load python/3.7-anaconda3

################################################################# Set Locations and Variables ##############################################################

#Set file locations
signal_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/merge_dir/master.signals.txt
command_dir=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/downsampling_variant_commands
output_dir=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/downsampling_variant_lists
stats_dir=/gpfs/alpine/med112/proj-shared/results/FOR_SUSIE

#Set populations
downsampled_pop=EUR
comparison_pop=AFR

#Set script locations
wrapper_script=/ccs/home/mconery/Pan-Trait_Colocalization/6_downsampling_analysis/make_downsampling_variant_lists_wrapper.py 
execution_script=/ccs/home/mconery/Pan-Trait_Colocalization/6_downsampling_analysis/execute_command_file.py

###########################################################################################################################################################

#Call python wrapper script with given parameters
python3 $wrapper_script -i $signal_file -c $command_dir -o $output_dir -s $stats_dir -n $downsampled_pop -t $comparison_pop

#Iterate over the command file and run the lines
ls $command_dir/*.downsampling_variant_list_commands.txt > $command_dir/make_downsampling_variant_commands_file_list.txt
declare -a cmd_array=()
while read line; do
	cmd_array+=("$line")
done < $command_dir/make_downsampling_variant_commands_file_list.txt

echo "JOB SIZE:"
echo "${#cmd_array[@]}"
for (( i=0; i<${#cmd_array[@]}; i++ )); do
	filename=$(basename ${cmd_array[$i]})
	filename=${filename%.downsampling_variant_list_commands.txt}
	echo jobname=$i for $filename
	srun -n1 -c1 --exclusive --job-name=$i python3 $execution_script ${cmd_array[$i]} &
done

wait
