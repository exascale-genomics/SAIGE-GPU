#!/bin/bash

#SBATCH -A med112
#SBATCH -N 60
#SBATCH -J EUR_make_matrix
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA//logs/count_filtered_AMR.out
#SBATCH -e /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA//logs/count_filtered_AMR.err

module load python/3.7-anaconda3

#################################################################### User set variables ###################################################################

ancestry=AMR

###########################################################################################################################################################

#Set script locations
execution_script=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/scripts/Pan-Trait_Colocalization.alex/2_matrix_making/execute_command.py
command_script=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/scripts/Pan-Trait_Colocalization.alex/4_signal_merging/count_filtered_snps_cmd.sh

#Set files and directories
command_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/count_filtered.$ancestry.cmds.txt 
susie_dir=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_"$ancestry"
out_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/filtered_out_dir/filtered_out.$ancestry.txt

#Call command making script
$command_script $susie_dir $command_file $out_file

#Iterate over the command file and run the lines
declare -a cmd_array=()
while read line; do
		cmd_array+=("$line")
done < $command_file

echo "JOB SIZE:"
echo "${#cmd_array[@]}"
counter=0
for (( i=0; i<${#cmd_array[@]}; i++ )); do
	srun -N1 -n1 -c4 --exclusive python3 $execution_script ${cmd_array[$i]} &
	counter=$((counter+1))
	if [ $counter -ge 480 ]; then # (32/4)*60
	    wait
        fi
done

wait
