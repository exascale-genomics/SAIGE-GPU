#!/bin/bash

#SBATCH -A med112
#SBATCH -N 4
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA//logs/finemapping.define_sig_loci.parallel.out
#SBATCH -e /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA//logs/finemapping.define_sig_loci.parallel.err
#SBATCH --distribution=cyclic:cyclic

module load python/3.7-anaconda3

################################################################# Set files and locations #################################################################

scripts_dir=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/scripts/Pan-Trait_Colocalization.alex/1_locus_defining #Directory corresponding to python code locations i.e. github repo
gwas_dir=/gpfs/alpine/med112/proj-shared/results/FOR_SUSIE #This is the new directory
locus_dir=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/defined_loci/ #Directory to hold locus outputs
final_file=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/final_output.grch37.high_thresh.broad.csv #Final output file (I recommend not putting it in the locus directory!)
primary_p_thresh=4.6e-11
secondary_p_thresh=5e-8
ancestry=META
align=grch37
exclusion_type=mhc_only

###########################################################################################################################################################

#### Run first fine-mapping excluding only mhc #####
chrom_bounds_file=$scripts_dir/"$align"_chr_bounds_"$exclusion_type".bed
declare -a cmd_array=()
#Loop over files in gwas_dir and execute code
#ls $gwas_dir/*META.GIA.KDI.txt.gz | while read file; do
for file in $gwas_dir/*META.GIA.KDI.txt.gz; do
	#Verify that file is for the right ancestry
	file=`basename $file`
	declare -a array=($(echo $file | tr "." " "))
	if [ ${array[1]} == $ancestry ]; then
		trait=${array[0]}
		out_file=$locus_dir/$align.$exclusion_type.$trait.$ancestry.sig_loci.csv
		#Check whether the individual trait already exists
		if [ ! -f "$out_file" ]; then
			echo $out_file
			#Send info to new shell script
			cmd="python3 $scripts_dir/define_sig_loci_for_trait.py -g $gwas_dir/$file -o $out_file -c $chrom_bounds_file -p $primary_p_thresh -s $secondary_p_thresh -n TRAIT.ANC.GIA.KDI.txt.gz"
			cmd_array+=("$cmd")
		fi
	fi
done

echo "JOB SIZE:"
echo "${#cmd_array[@]}"
counter=0
for i in "${cmd_array[@]}"; do
	echo $i
	srun -N1 -n1 -c2 --exclusive $i &
	counter=$((counter+1))
        if [ $counter -ge 256 ]; then # (32/4)*60
            wait
	    counter=0
        fi
done

wait

#Submit final job contingent on other jobs finishing
srun -n1 -c1 --exclusive python3 $scripts_dir/combine_sig_loci.py -l $locus_dir -o $final_file -a META -f $align.$exclusion_type.TRAIT.ANC.sig_loci.csv
