#!/bin/bash

#SBATCH -A med112
#SBATCH -N 50
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA//logs/finemapping.fine-map.EUR.parallel.run1_sorted.Anurag.out
#SBATCH -e /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA//logs/finemapping.fine-map.EUR.parallel.run1_sorted.Anurag.err
#SBATCH --distribution=cyclic:cyclic

module load python/3.7-anaconda3
module load r/4.0.3-py3



############################################################### Set run-specific variables ################################################################

ancestry=EAS #Also change in output file names
user=Anurag
batch=small

###########################################################################################################################################################

#Script locations
wrapper_script=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/scripts/Pan-Trait_Colocalization.alex/3_finemap/check_and_map_loci.py 
execution_script=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/scripts/Pan-Trait_Colocalization.alex/2_matrix_making/execute_command.py

#Files and directories
locus_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/final_output.grch37.high_thresh.broad."$batch"_subset.csv
matrix_dir=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/matrix_dir_"$ancestry" #Change to desired matrix directory that holds matrices and map files
susie_dir=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_"$ancestry" #Output directory
stats_dir=/gpfs/alpine/med112/proj-shared/results/FOR_SUSIE #Directory of all summary stats
command_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/fine-map.$ancestry.$user.cmds.txt #Change user for file
stats_name_conv="TRAIT.ANC.GIA.KDI.txt.gz"

#Call python wrapper script with given parameters
python3 $wrapper_script -l $locus_file -c $command_file -o $susie_dir -m $matrix_dir -s $stats_dir -a $ancestry -n $stats_name_conv

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
