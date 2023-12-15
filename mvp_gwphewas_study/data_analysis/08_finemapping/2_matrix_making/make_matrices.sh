#!/bin/bash

#SBATCH -A med112
#SBATCH -N 60
#SBATCH -J EUR_make_matrix
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA//logs/finemapping.make_matrices.parallel.EUR.new.run1_sorted.out
#SBATCH -e /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA//logs/finemapping.make_matrices.parallel.EUR.new.run1_sorted.err
#SBATCH --distribution=cyclic:cyclic

module load python/3.7-anaconda3

################################################################# Set files and locations #################################################################

#Set ancestry (Will need to run for each pop!)
ancestry=EUR 

#Script locations
wrapper_script=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/scripts/Pan-Trait_Colocalization.alex/2_matrix_making/make_matrices.py 
execution_script=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/scripts/Pan-Trait_Colocalization.alex/2_matrix_making/execute_command.py

#Files and directories
plink_dir=/gpfs/alpine/med112/proj-shared/task0101113/output/pheCodes/inputs/bgen #Actual plink file directory containing the pgen files
matrix_dir=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/matrix_dir_"$ancestry"
individual_dir=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/phenotype_individual_dir
locus_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/final_output.grch37.high_thresh.broad.csv #This list is ascending size-sorted
registry_file=/gpfs/alpine/med112/proj-shared/data/genetic_annotations/MVP_genetic_data_chunks.txt #Registry of plink file chunks
command_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/make_matrices.cmds.txt #Fill in location for temporary file of commands
ancestry_file_format="PHENO.ANC.individuals.txt"
max_dist=1000000000

###########################################################################################################################################################

#Call python wrapper script with given parameters
python3 $wrapper_script -l $locus_file -o $matrix_dir -a $ancestry -i $individual_dir -r $registry_file -c $command_file -p $plink_dir -f $ancestry_file_format -d $max_dist

#Iterate over the command file and run the lines
declare -a cmd_array=()
while read line; do
		cmd_array+=("$line")
done < $command_file

echo "JOB SIZE:"
echo "${#cmd_array[@]}"
for (( i=0; i<${#cmd_array[@]}; i++ )); do
	filename=$(basename ${cmd_array[$i]})
	filename=${filename%.ped}
	echo jobname=$i for $filename
	srun -n1 -c32 --exclusive --job-name=$i python3 $execution_script ${cmd_array[$i]} &
	if [ $i -ge 480 ]; then # (32/4)*60
	    wait
        fi
done

wait
