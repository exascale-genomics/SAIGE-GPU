#!/bin/bash
#SBATCH -A bif154
#SBATCH -J saige_step1_cpu
#SBATCH -o run_step1_saige_all_traits.cpu.log
#SBATCH -N 10
#SBATCH -p batch
#SBATCH -t 1:00:00
#SBATCH -C nvme

##################################################################################################################
################################################### Set Needed Variables #########################################
##################################################################################################################
#Input directories
plink_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes"
pheno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes"
out_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/CPU_test"

#Set up step1 single-trait script location
saige_step1_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/run_step1_saige_single_trait.cpu.sh"

#Set population names
populations=("AFR")
trait_types=("binary")

##################################################################################################################
#################################################### Run SAIGE Step 1 ############################################
##################################################################################################################

#Get total number of nodes
NNODES=${SLURM_NNODES}
#Get maximum number of parallel jobs
max_parallel=$NNODES

cmds=() #Create empty commands vector
#Loop over the populations, create the mapping commands for all phenotypes, and launch them
for pop in ${populations[@]}; do
	plink_file_pre=$plink_dir/$pop/merged.plink.maf.ld
	plink_file_post=/mnt/bb/${USER}/$pop
	#sbcast genotype files to local scratch
	sbcast -pf $plink_file_pre.bim $plink_file_post.bim
	sbcast -pf $plink_file_pre.fam $plink_file_post.fam
	sbcast -pf $plink_file_pre.bed $plink_file_post.bed
	#Change into output directory
	pop_out_dir=$out_dir/$pop
	mkdir -p $pop_out_dir
	cd $pop_out_dir
	touch $pop_out_dir/timing.csv  # Create timing file
	#Create commands to process step 1 for all phenotypes of the population
	#Loop over the trait types
	for trait_type in ${trait_types[@]}; do
		#Make a list of the phenotype files for the given trait type and population
		ls $pheno_dir/$pop/Pheno*.$trait_type.txt > $pheno_dir/$pop/$trait_type.phenotype_files.txt
		#Loop over the files and create the commands
		while read file; do
			temp=$(basename $file)
			phenotype=${temp%.$trait_type.txt}
			#Check if output files exist and add command if they do not
			cmds+=("$saige_step1_script $plink_file_post $file $phenotype $pop_out_dir $trait_type")
		done < $pheno_dir/$pop/$trait_type.phenotype_files.txt
	done
done

#Run all the commands for the populations
job_size=${#cmds[@]}
for (( i=0; i<10; i++ )); do
	j=$(( $i + 1 ))
	srun -N1 -n1 -c56 --exclusive -C nvme --job-name=$j ${cmds[$i]} & 
	if [ $(($j % $max_parallel)) == 0 ] || [ $j == $job_size ]; then 
	    	wait
        fi
done