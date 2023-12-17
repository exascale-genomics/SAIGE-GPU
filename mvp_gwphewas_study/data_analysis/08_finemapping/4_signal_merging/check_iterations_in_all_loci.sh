#!/bin/bash

#SBATCH -A med112
#SBATCH -N 50
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/logs/finemapping.iter_check.AFR.parallel.run1_sorted.Alex.out
#SBATCH -e /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/logs/finemapping.iter_check.AFR.parallel.run1_sorted.Alex.err
#SBATCH --distribution=cyclic:cyclic

module load python/3.7-anaconda3
module load r/4.0.3-py3


############################################################### Set run-specific variables ################################################################

user=Alex
ancestry=AFR

###########################################################################################################################################################

#Script locations
wrapper_script=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/scripts/Pan-Trait_Colocalization.alex/4_signal_merging/check_iterations_in_all_loci.py
execution_script=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/scripts/Pan-Trait_Colocalization.alex/2_matrix_making/execute_command.py

#Files and directories
ANC_dir=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_"$ancestry"
command_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/check_iter.$user.$ancestry.cmds.txt #Change user for file

#Call python wrapper script with given parameters
python3 $wrapper_script -m $command_file -a $ANC_dir

#Iterate over the command file and run the lines
declare -a cmd_array=()
while read line; do
		cmd_array+=("$line")
done < $command_file

echo "JOB SIZE:"
echo "${#cmd_array[@]}"
counter=0
for (( i=0; i<${#cmd_array[@]}; i++ )); do
	filename=${cmd_array[$i]}
        filename=${filename%___-o___*}
	filename=${filename#*___-p___}
	echo jobname=$i for $filename
	srun -n1 -c32 --exclusive --job-name=$i python3 $execution_script ${cmd_array[$i]} &
        counter=$((counter+1))
        if [ $counter -ge 200 ]; then # (32/4)*60
            wait
            counter=0
        fi
done

wait
