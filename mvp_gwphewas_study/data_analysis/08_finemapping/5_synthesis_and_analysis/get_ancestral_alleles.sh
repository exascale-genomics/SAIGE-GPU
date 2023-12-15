###############################Define Directories###################################

#Set directories and locations
work_dir=/project/voight_GWAS/mconery/downstream_analyses
vcf_dir=/project/voight_datasets_01/1kg/phaseIII_2013
int_dir=$work_dir/high_pip_filtered_vcfs

####################################################################################

#Load modules
module load bcftools/1.15.1

#Loop over chromosomes
for chromo in {1..22}; do 
	bsub -R "rusage[mem=24G]" -M 25G "bcftools view -G -T $work_dir/high_pip_targets.for_vcf.txt $vcf_dir/ALL.chr$chromo.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep ^$chromo | sed 's/AC=[^$]*;AA/AA/g' | sed 's/|||[^$]*//g' > $int_dir/chr$chromo.ancestral.txt"
done



