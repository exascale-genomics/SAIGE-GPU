#!/bin/bash

#SBATCH -A med112
#SBATCH -N 60
#SBATCH -J EUR_make_matrix
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA//logs/matrix_qc_AMR.out
#SBATCH -e /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA//logs/matrix_qc_AMR.err

module load python/3.7-anaconda3

################################################################# Set files and locations #################################################################

#Script locations
#Files and directories
execution_script=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/scripts/Pan-Trait_Colocalization.alex/2_matrix_making/execute_command.py
command_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/make_qc_matrices.AMR.cmds.txt #Fill in location for temporary file of commands

#Usage parameters
#max_ram=25 #Currently not being used

###########################################################################################################################################################

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
