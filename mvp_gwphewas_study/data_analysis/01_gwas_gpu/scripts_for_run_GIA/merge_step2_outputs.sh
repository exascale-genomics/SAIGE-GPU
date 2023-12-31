#!/usr/bin/bash
# merge all results from step2 downsample run
work_dir=.//output/GIA_ANC_Run
phe_dir=./results/PostGWAS_Analysis/finemaping/GIA/downsampling_variant_lists
group="EUR"
output_path="$work_dir/downsample_summary_stats"
header=`head -n1 ./output/GIA_ANC_Run/Weight_Max_INT/EUR/step2/chr12.10.dose.bgen.Weight_Max_INT.EUR.downsample_orig.assoc.gpu.txt`

# find all downsampled files
for i in $phe_dir/*.EUR.downsampling_variants.txt; do
	pheno=`basename $i .EUR.downsampling_variants.txt`
	echo $pheno
	step2_path="$work_dir/$pheno/$group/step2"
	summary_out=$output_path/$pheno.$group.downsample.summary.txt
	echo $header > $summary_out
	cat $step2_path/*.downsample.assoc.gpu.txt | grep -v "CHR" >> $summary_out;
done
