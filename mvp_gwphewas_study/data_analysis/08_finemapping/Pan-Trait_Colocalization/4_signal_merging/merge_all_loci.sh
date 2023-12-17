#!/bin/bash

#SBATCH -A med112
#SBATCH -N 50
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/logs/finemapping.merge.parallel.run1_sorted.Alex.out
#SBATCH -e /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/logs/finemapping.merge.parallel.run1_sorted.Alex.err
#SBATCH --distribution=cyclic:cyclic

source /ccs/home/arodriguez/.bashrc
conda activate R_4.0_mkl


############################################################### Set run-specific variables ################################################################

user=Alex

###########################################################################################################################################################

#Script locations
wrapper_script=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/scripts/Pan-Trait_Colocalization.alex/4_signal_merging/merge_all_loci.py
execution_script=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/scripts/Pan-Trait_Colocalization.alex/2_matrix_making/execute_command.py

#Files and directories
locus_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/final_output.grch37.high_thresh.broad.csv
AFR_dir=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_AFR
AMR_dir=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_AMR
EAS_dir=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_EAS
EUR_dir=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_EUR
output_dir=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/merge_dir
command_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/merge.$user.cmds.txt #Change user for file
master_file=$output_dir/master.signals.txt
variant_file=$output_dir/master.variants.txt

#Set variables
susie_name_conv="TRAIT.CHR.START.END.ANC.rds"
sig_thresh=5e-8
sig_type=gwas
purity_thresh=0.1
cut_height=0.9

#Call python wrapper script with given parameters
python3 $wrapper_script -l $locus_file -m $command_file -o $output_dir -a $AFR_dir -b $AMR_dir -c $EAS_dir -d $EUR_dir -t $sig_type -s $sig_thresh -p $purity_thresh -u $cut_height -n $susie_name_conv -e True

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
	srun -n1 -N1 -c2 --exclusive --job-name=$i python3 $execution_script ${cmd_array[$i]} &
        counter=$((counter+1))
        if [ $counter -ge 1000 ]; then # (32/4)*60
            wait
            counter=0
        fi
done

wait

#Merge signal-level outputs into a single file
>$master_file
rm $master_file
for file in $output_dir/*.MERGE.txt; do 
	if test -f "$master_file"; then
    		awk 'NR>1' $file >> $master_file
	else
		cat $file > $master_file
	fi
done

#Merge variant-level outputs into a single file
>$variant_file
rm $variant_file
for file in $output_dir/*.MERGE_VARIANT.txt; do 
	if test -f "$variant_file"; then
    		awk 'NR>1' $file >> $variant_file
	else
		cat $file > $variant_file
	fi
done