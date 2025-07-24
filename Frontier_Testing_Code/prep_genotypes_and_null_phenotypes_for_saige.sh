#!/bin/bash
#SBATCH -A bif154
#SBATCH -J prep_geno
#SBATCH -o prep_genotypes_and_null_phenotypes_for_saige.log
#SBATCH -N 8
#SBATCH -p batch
#SBATCH -t 1:00:00

#Load plink tools and conda environment
source ~/.bashrc
conda activate simulate_pheno

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################
#Set directories
inp_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/simulated_genotypes"
geno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes"
pheno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes"

#Set script locations
null_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/simulate_null_phenotypes.py"
simulate_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/simulate_phenotype.py"
decode_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/3_Fine-Mapping_Robustness_Testing/execute_command.py"

#Get number of nodes
num_nodes="$SLURM_NNODES"

################################################################################################################
############################## Filter for MAF and Make LD Pruning Filters ######################################
################################################################################################################
##Create empty command array
#cmds=()
##Loop over two populations
#for pop in AFR EUR; do 
#	#Change into population output directory
#	mkdir -p $geno_dir/$pop
#	cd $geno_dir/$pop
#	#Get keep file location
#	keep_file=$inp_dir/$pop/$pop.400K_random_sample.tsv
#	#Loop the chromosomes
#	for chromo in {1..22}; do
#		#Set plink file path
#		plink_inp=$inp_dir/$pop/chr"$chromo".plink.tag
#		plink_out=$geno_dir/$pop/chr"$chromo".plink.maf
#		if [ $pop == "EUR" ]; then 
#			cmds+=("plink2 --bfile $plink_inp --keep $keep_file --set-all-var-ids "@:#:\$r:\$a" --maf 0.01 --indep-pairwise 1000 80 0.1 --make-bed --out $plink_out")
#		else
#			cmds+=("plink2 --bfile $plink_inp --keep $keep_file --set-all-var-ids "@:#:\$r:\$a" --maf 0.01 --indep-pairwise 1000 80 0.05 --make-bed --out $plink_out")
#		fi
#	done
#done
#
##Launch commands
#job_size=${#cmds[@]}
#echo "JOB SIZE: $job_size"
#for (( i=0; i<${#cmds[@]}; i++ )); do
#	j=$(( $i + 1 ))
#	srun -n1 -N1 --ntasks-per-node=1 --job-name=$j ${cmds[$i]} &
#	if [ $(($j % $num_nodes)) == 0 ] || [ $j == $job_size ]; then 
#	    wait
#        fi
#done
#wait
#
################################################################################################################
################################################# Prune the Variants ###########################################
################################################################################################################
##Create empty command array
#cmds=()
##Loop over two populations
#for pop in AFR EUR; do 
#	#Change into population output directory
#	mkdir -p $geno_dir/$pop
#	cd $geno_dir/$pop
#	#Get keep file location
#	keep_file=$inp_dir/$pop/$pop.400K_random_sample.tsv
#	#Loop the chromosomes
#	for chromo in {1..22}; do
#		#Set plink file path
#		plink_inp=$inp_dir/$pop/chr"$chromo".plink.tag
#		plink_out=$geno_dir/$pop/chr"$chromo".plink.maf
#		#Launch plink to prune variants
#		cmds+=("plink2 --bfile $plink_out --extract $plink_out.prune.in --make-bed --out $plink_out.ld")
#	done
#done
#
##Launch commands
#job_size=${#cmds[@]}
#echo "JOB SIZE: $job_size"
#for (( i=0; i<${#cmds[@]}; i++ )); do
#	j=$(( $i + 1 ))
#	srun -n1 -N1 --ntasks-per-node=1 --job-name=$j ${cmds[$i]} &
#	if [ $(($j % $num_nodes)) == 0 ] || [ $j == $job_size ]; then 
#	    wait
#        fi
#done
#wait
#
################################################################################################################
################################################# Make Merged Files ############################################
################################################################################################################
##Create empty command array
#cmds=()
##Loop over the populations
#for pop in AFR EUR; do
#	cd $geno_dir/$pop
#	#Create merge list and merge files
#	ls $geno_dir/$pop/chr*.plink.maf.ld.bim | sed 's/.bim//g' > $geno_dir/$pop/$pop.mergelist.txt
#	cmds+=("plink2 --allow-no-sex --make-bed --pmerge-list $geno_dir/$pop/$pop.mergelist.txt bfile --out $geno_dir/$pop/merged.plink.maf.ld")
#done
#
##Launch commands
#job_size=${#cmds[@]}
#echo "JOB SIZE: $job_size"
#for (( i=0; i<${#cmds[@]}; i++ )); do
#	j=$(( $i + 1 ))
#	srun -n1 -N1 --ntasks-per-node=1 --job-name=$j ${cmds[$i]} &
#	if [ $(($j % $num_nodes)) == 0 ] || [ $j == $job_size ]; then 
#	    wait
#        fi
#done
#wait

################################################################################################################
################################################ Calculate PCs #################################################
################################################################################################################
##Create empty command array
#cmds=()
##Loop over the populations
#for pop in AFR EUR; do
#	mkdir -p $pheno_dir/$pop
#	cd $pheno_dir/$pop
#	#Calculate PCs
#	cmds+=("plink2 --bfile $geno_dir/$pop/merged.plink.maf.ld --pca approx 10 --out $geno_dir/$pop/merged.plink.maf.ld.pcs")
#done
#
##Launch commands
#job_size=${#cmds[@]}
#echo "JOB SIZE: $job_size"
#for (( i=0; i<${#cmds[@]}; i++ )); do
#	j=$(( $i + 1 ))
#	srun -n1 -N1 --ntasks-per-node=1 --job-name=$j ${cmds[$i]} &
#	if [ $(($j % $num_nodes)) == 0 ] || [ $j == $job_size ]; then 
#	    wait
#        fi
#done
#wait
#
################################################################################################################
############################################ Simulate Null Phenotypes ##########################################
################################################################################################################
##Create empty command array
#cmds=()
##Loop over the populations
#for pop in AFR EUR; do
#	cd $pheno_dir/$pop
#	cmds+=("python $null_script -i $geno_dir/$pop/merged.plink.maf.ld.pcs.eigenvec -o $pheno_dir/$pop/merged.maf.ld.continuous.txt -n 1000 -t continuous")
#	cmds+=("python $null_script -i $geno_dir/$pop/merged.plink.maf.ld.pcs.eigenvec -o $pheno_dir/$pop/merged.maf.ld.binary.txt -n 1000 -t binary")
#done
#
##Launch commands
#job_size=${#cmds[@]}
#echo "JOB SIZE: $job_size"
#for (( i=0; i<${#cmds[@]}; i++ )); do
#	j=$(( $i + 1 ))
#	srun -n1 -N1 --ntasks-per-node=1 --job-name=$j ${cmds[$i]} &
#	if [ $(($j % $num_nodes)) == 0 ] || [ $j == $job_size ]; then 
#	    wait
#        fi
#done
#wait

################################################################################################################
################################################# Split Phenotypes #############################################
################################################################################################################
##Make commands to split the files
#cmds=()
##Loop over the populations
#for pop in AFR EUR; do for trait_type in binary continuous; do cd $pheno_dir/$pop; for i in {1..1000}; do
#		pheno_col=$(( $i + 12))
#		if [ ! -f "$pheno_dir"/"$pop"/Pheno"$i"."$trait_type".txt ]; then 
#			cmds+=("python $decode_script -c cut___-f___1-12,"$pheno_col"___"$pheno_dir"/"$pop"/merged.maf.ld."$trait_type".txt___>___"$pheno_dir"/"$pop"/Pheno"$i"."$trait_type".txt")
#		fi
#done; done; done
#
##Calculate number of possible jobs
#job_size=${#cmds[@]}
#num_parallel_jobs=$((($SLURM_CPUS_ON_NODE / 4) * num_nodes))
#echo "JOB SIZE: $job_size"
##Launch commands
#for (( i=0; i<${#cmds[@]}; i++ )); do
#	j=$(( $i + 1 ))
#	srun -n1 -c1 -N1 --exclusive --job-name=$j ${cmds[$i]} &
#	if [ $(($j % $num_parallel_jobs)) == 0 ] || [ $j == $job_size ]; then 
#	    wait
#        fi
#done
#wait

################################################################################################################
################################## Replace the Phenotype Name in the Split File ################################
################################################################################################################
##Make commands to split the files
#cmds=()
##Loop over the populations
#for pop in AFR EUR; do for trait_type in binary continuous; do cd $pheno_dir/$pop; for i in {1..1000}; do
#			cmds+=("python $decode_script -c sed___-i___'s/Pheno$i/PHENO/g'___"$pheno_dir"/"$pop"/Pheno"$i"."$trait_type".txt")
#done; done; done
#
##Calculate number of possible jobs
#job_size=${#cmds[@]}
#num_parallel_jobs=$((($SLURM_CPUS_ON_NODE / 4) * num_nodes))
#echo "JOB SIZE: $job_size"
##Launch commands
#for (( i=0; i<${#cmds[@]}; i++ )); do
#	j=$(( $i + 1 ))
#	srun -n1 -c1 -N1 --exclusive --job-name=$j ${cmds[$i]} &
#	if [ $(($j % $num_parallel_jobs)) == 0 ] || [ $j == $job_size ]; then 
#	    wait
#        fi
#done
#wait

################################################################################################################
############################################## Make Step-2 BGEN Files ##########################################
################################################################################################################
##Create empty command array
#cmds=()
##Loop over two populations
#for pop in EUR; do
#  #Change into population output directory
#  cd $geno_dir/$pop
#  #Loop the chromosomes
#  for chromo in {1..22}; do
#    #Set input plink file prefix (without extension)
#    plink_inp=$geno_dir/$pop/chr"$chromo".plink.maf
#    #Set output bgen prefix
#    bgen_out=$geno_dir/$pop/chr"$chromo".plink.maf.bgen
#    #Command to convert to BGEN v1.2 with 8-bit compression
#    cmds+=("plink2 --bfile $plink_inp --export bgen-1.2 bits=8 --out $plink_inp")
#  done
#done
#
##Launch commands
#job_size=${#cmds[@]}
#echo "JOB SIZE: $job_size"
#for (( i=0; i<${#cmds[@]}; i++ )); do
#	j=$(( $i + 1 ))
#	srun -n1 -N1 --ntasks-per-node=1 --job-name=$j ${cmds[$i]} &
#	if [ $(($j % $num_nodes)) == 0 ] || [ $j == $job_size ]; then 
#	    wait
#       fi
#done
#wait

################################################################################################################
######################################### Make BGEN Index Files (.bgen.bgi) ####################################
################################################################################################################
#cmds=()
##Loop over two populations
#for pop in EUR; do
#  cd $geno_dir/$pop
#  #Loop the chromosomes
#  for chromo in {1..22}; do
#    bgen_file=$geno_dir/$pop/chr"$chromo".plink.maf.bgen
#    #Command to create index for each BGEN file
#    cmds+=("bgenix -g $bgen_file -index")
#  done
#done
#
##Launch index commands
#job_size=${#cmds[@]}
#echo "JOB SIZE: $job_size"
#for (( i=0; i<${#cmds[@]}; i++ )); do
#  j=$(( $i + 1 ))
#  srun -n1 -N1 --ntasks-per-node=1 --job-name=$j ${cmds[$i]} &
#  if [ $(($j % $num_nodes)) == 0 ] || [ $j == $job_size ]; then
#    wait
#  fi
#done
#wait

################################################################################################################
############################ Make BGEN Variant Lists (for Multi-Threading SAIGE Step 2) ########################
################################################################################################################
cmds=()
#Loop over two populations
for pop in AFR EUR; do
  cd $geno_dir/$pop
  #Loop the chromosomes
  for chromo in {1..22}; do
    bim_file=$geno_dir/$pop/chr"$chromo".plink.maf.bim
    var_file=$geno_dir/$pop/chr"$chromo".plink.maf.var_list.txt
    #Command to create index for each BIM file
    awk '{print $2}' $bim_file > $var_file
  done
done

#Launch index commands
job_size=${#cmds[@]}
echo "JOB SIZE: $job_size"
for (( i=0; i<${#cmds[@]}; i++ )); do
  j=$(( $i + 1 ))
  srun -n1 -N1 --ntasks-per-node=1 --job-name=$j ${cmds[$i]} &
  if [ $(($j % $num_nodes)) == 0 ] || [ $j == $job_size ]; then
    wait
  fi
done
wait

################################################################################################################
############################################### Make Sample Lists ##############################################
################################################################################################################
##Make the sample lists needed for step 2
#for pop in AFR EUR; do
#	awk 'NR > 1{print $2}' $geno_dir/$pop/merged.plink.maf.ld.psam > $geno_dir/$pop/sample_list.txt
#done