#!/usr/bin/bash
samples=./results/DOWNSAMPLING/Samples.EUR.AFR_downsample.txt
phe_dir=./results/PostGWAS_Analysis/finemaping/GIA/downsampling_variant_lists
work_dir=./output/GIA_ANC_Run/
rscript_file=./batch/pre-gwas/gwPheWAS-Summit/scripts_for_run_GIA/create_step2_rds_downsample.r
group="EUR"

module load r/4.0.5

# create the rds file for each phenotype to be submitted
# only create jobs for the chromosomes listed in the downsampling_variants files
for i in $phe_dir/*.EUR.downsampling_variants.txt; do
        #echo $i
        pheno=`basename $i .EUR.downsampling_variants.txt`
	chrms=`cut -f1 -d: $i | sort | uniq | paste -sd, -`
	outrds=$work_dir/$pheno/$group/submit/step2/submit.step2.downsampled.$pheno.$group.rds
	echo "Rscript $rscript_file $pheno $group $outrds \"$chrms\" $i" >> command_list.txt
	#Rscript $rscript_file $pheno $group $outrds \"$chrms\" $i
done

