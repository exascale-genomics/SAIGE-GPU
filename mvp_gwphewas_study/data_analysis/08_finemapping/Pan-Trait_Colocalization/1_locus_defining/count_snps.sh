#!/bin/bash

#SBATCH -A med112
#SBATCH -N 20
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA//logs/finemapping.count_snps.parallel.EUR.run2.out
#SBATCH -e /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA//logs/finemapping.count_snps.parallel.EUR.run2.err

module load python/3.7-anaconda3

################################################################# Set files and locations #################################################################

#Script locations
wrapper_script=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/scripts/Pan-Trait_Colocalization.alex/1_locus_defining/count_snps.py 
execution_script=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/scripts/Pan-Trait_Colocalization.alex/2_matrix_making/execute_command.py

#Files and directories
plink_dir=/gpfs/alpine/med112/proj-shared/task0101113/output/pheCodes/inputs/bgen 
matrix_dir=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/count_snps_dir #Change to desired matrix directory that holds matrices and map files
locus_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/final_output.test.csv 
#locus_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/final_output.subset.csv
registry_file=/gpfs/alpine/med112/proj-shared/data/genetic_annotations/MVP_genetic_data_chunks.txt
ancestry_file=/gpfs/alpine/med112/proj-shared/GIA/output/projections/smartpca_loadings_projection/release4.mvp.gia.EUR.txt #Will need to run for each pop!
command_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/count_snps.cmds.txt 
ancestry_file_format="XXX.XXX.XXX.ANC.txt"

#Usage parameters
#max_ram=25 #Currently not being used

###########################################################################################################################################################

#Call python wrapper script with given parameters
python3 $wrapper_script -l $locus_file -o $matrix_dir -a $ancestry_file -r $registry_file -c $command_file -p $plink_dir -f $ancestry_file_format

#Iterate over the command file and run the lines
declare -a cmd_array=()
while read line; do
		cmd_array+=("$line")
done < $command_file

echo "JOB SIZE:"
echo "${#cmd_array[@]}"
for i in "${cmd_array[@]}"; do
	srun -n1 -c32 --exclusive python3 $execution_script $i &
done

wait

#Count SNPs in each pvar file and print to console
wc -l $matrix_dir/*.pvar
