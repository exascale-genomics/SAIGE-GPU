# this Bash script will create the downsampled pheno files for the selected phecodes in the $phe_dir 

samples=./results/DOWNSAMPLING/Samples.EUR.AFR_downsample.txt
phe_dir=./results/DOWNSAMPLING/SNPLIST_dwnsmpl_AFR
work_dir=./output/GIA_ANC_Run/
group="EUR"

for i in $phe_dir/*.EUR.downsampling_variants.txt; do
	echo $i
	pheno=`basename $i .EUR.downsampling_variants.txt`
	ophenotype_file=$work_dir/$pheno/$group/inputs/PhenoFile.$pheno.$group.txt
	pheno_output=$work_dir/$pheno/$group/inputs/downsampled.PhenoFile.$pheno.$group.txt
	head -n1 $ophenotype_file > $pheno_output
	cut -f1 $samples | grep -f - $ophenotype_file >> $pheno_output 
done
