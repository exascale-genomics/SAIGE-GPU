#!/usr/bin/bash
work_dir=.//output/GIA_ANC_Run
phe_dir=./results/DOWNSAMPLING/SNPLIST_dwnsmpl_AFR
group="EUR"
submit_file=$work_dir/SUBMIT/downsampled.submit.step1.parallel.sh
step1_path=./tools/saige_20220326/SAIGE-DOE/extdata/step1_fitNULLGLMM.R
plink_files=./output/pheCodes/inputs/genotypes/re-run/20200917.GenotypeData.Release4.mvpcoreid.ld.maf.05
quant_file=./data/phenotype_data/processed/GIA/gia_quant.txt
binary_file=./data/phenotype_data/processed/GIA/gia_binary.txt
female_binary_codes=./data/phenotype_data/processed/PhenotypeList.Binary.Female.Traits.v2.txt
male_binary_codes=./data/phenotype_data/processed/PhenotypeList.Binary.Male.Traits.v2.txt
female_quant_codes=./data/phenotype_data/processed/PhenotypeList.Quantitative.Female.Traits.txt

echo "#!/bin/bash" > $submit_file
echo "#BSUB -nnodes 1" >> $submit_file
echo "#BSUB -W 4:00" >> $submit_file
echo "#BSUB -q batch-hm" >> $submit_file
echo "#BSUB -P MED112" >> $submit_file
echo "#BSUB -o all.s1.EUR.downsampled.gpu.stdout" >> $submit_file
echo "#BSUB -e all.s1.EUR.downsampled.gpu.stderr" >> $submit_file
echo "#BSUB -J all.s1.down.EUR" >> $submit_file
echo "#BSUB -alloc_flags nvme" >> $submit_file
echo "" >> $submit_file
echo "module load python/2.7.15-anaconda2-5.3.0" >> $submit_file
echo "module load r/4.0.5" >> $submit_file
echo "module load gcc/9.1.0" >> $submit_file
echo "module load cmake" >> $submit_file
echo "module load openblas" >> $submit_file
echo "module load cuda/11.0.2" >> $submit_file
echo "module load spectrum-mpi/10.4.0.3-20210112" >> $submit_file
echo "" >> $submit_file

# find all downsampled files
for i in $phe_dir/*.EUR.downsampling_variants.txt; do
	echo $i
	pheno=`basename $i .EUR.downsampling_variants.txt`
	pheno_file=$work_dir/$pheno/$group/inputs/downsampled.PhenoFile.$pheno.$group.txt
	output_prefix=$work_dir/$pheno/$group/step1/downsampled.phewas.ld.maf.$pheno.$group.out.gpu
	log_path=$work_dir/$pheno/$group/step1/downsampled.phewas.ld.maf.$pheno.$group.out.gpu.log
	traitType=binary
	if [ "$(grep -c "^$pheno$" $quant_file)" -ge 1 ]; then
		traitType=quantitative
	fi

	female_params="--FemaleOnly TRUE --FemaleCode 2 --sexCol Gender "
	male_params="--MaleOnly TRUE --MaleCode 1 --sexCol Gender "
	gender_params=""
	if [ "$(grep -c "^$pheno$" $female_binary_codes)" -ge 1 ]; then
                gender_params=$female_params
        fi
	if [ "$(grep -c "^$pheno$" $female_quant_codes)" -ge 1 ]; then
                gender_params=$female_params
        fi
        if [ "$(grep -c "^$pheno$" $male_binary_codes)" -ge 1 ]; then
                gender_params=$male_params
        fi


	# use 2 GPUs for this size of samples (similar to AFR on batch-hm)
	echo "jsrun -n 2 -a1 -c1 -g1 Rscript $step1_path --plinkFile=$plink_files --phenoFile=$pheno_file --outputPrefix=$output_prefix --phenoCol=Phenotype --covarColList=Age,Gender,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --qCovarColList=Gender --sampleIDColinphenoFile=MVPCore_ID --traitType=$traitType --nThreads=1 --minMAFforGRM 0.05 --maxiterPCG 500 --maxiter 25 --LOCO TRUE $gender_params > $log_path &" >> $submit_file

done
echo "wait" >> $submit_file
