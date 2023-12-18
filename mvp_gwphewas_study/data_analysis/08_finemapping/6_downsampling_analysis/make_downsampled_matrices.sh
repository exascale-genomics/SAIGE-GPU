#!/bin/bash

#SBATCH -A med112
#SBATCH -N 1
#SBATCH -J EUR_make_matrix
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/logs/finemapping.make_downsampled_matrices.parallel.out
#SBATCH -e /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/logs/finemapping.make_downsampled_matrices.parallel.err
#SBATCH --distribution=cyclic:cyclic

module load python/3.7-anaconda3

################################################################# Set files and locations #################################################################

#Set Global Variables
down_anc=EUR
ref_anc=AFR
max_dist=1000000000

#Script locations
locus_script=/ccs/home/mconery/Pan-Trait_Colocalization/6_downsampling_analysis/define_downsampling_loci.py
wrapper_script=/ccs/home/mconery/Pan-Trait_Colocalization/2_matrix_making/make_matrices.py 
execution_script=/ccs/home/mconery/Pan-Trait_Colocalization/2_matrix_making/execute_command.py

#Files and directories
out_dir=/gpfs/alpine/proj-shared/med112/results/DOWNSAMPLING
locus_file=$out_dir/downsampling_loci.grch37.high_thresh.broad.csv #This list is ascending size-sorted
trait_file=$out_dir/downsampled_traits.tsv
signal_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/merge_dir/master.signals.txt
plink_dir=/gpfs/alpine/med112/proj-shared/task0101113/output/pheCodes/inputs/bgen #Actual plink file directory containing the pgen files
matrix_dir=$out_dir/matrix_dir_"$down_anc"_to_"$ref_anc"
individual_dir=$out_dir/phenotype_individuals_"$down_anc"_to_"$ref_anc"
registry_file=/gpfs/alpine/med112/proj-shared/data/genetic_annotations/MVP_genetic_data_chunks.txt #Registry of plink file chunks
command_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/make_downsampled_matrices."$down_anc"_to_"$ref_anc".cmds.txt #Fill in location for temporary file of commands
ancestry_file_format="PHENO.ANC.individuals.txt"

###########################################################################################################################################################

##Call locus defining script
#python3 $locus_script -i $signal_file -o $locus_file -p $trait_file -n $down_anc -t $ref_anc

#Make needed directories
mkdir -p $individual_dir
mkdir -p $matrix_dir

##Make lists of variants
#declare -a trait_cmd_array=()
#while read trait; do
#	gwas_trait_file=/gpfs/alpine/proj-shared/med112/task0101113/output/GIA_ANC_Run/$trait/EUR/inputs/downsampled.PhenoFile.$trait.EUR.txt
#	trait_cmd_array+=("awk 'NR > 1 {print \$1"\"\\t"\"\$1}' $gwas_trait_file > $individual_dir/$trait.EUR.individuals.txt")
#done < $trait_file
#echo "JOB SIZE:"
#echo "${#trait_cmd_array[@]}"
#counter=0
#for (( i=0; i<${#trait_cmd_array[@]}; i++ )); do
#	eval ${trait_cmd_array[$i]}
#done

#Call python wrapper script with given parameters
python3 $wrapper_script -l $locus_file -o $matrix_dir -a $down_anc -i $individual_dir -r $registry_file -c $command_file -p $plink_dir -f $ancestry_file_format -d $max_dist

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
