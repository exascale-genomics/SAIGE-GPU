# create the SNP list 
snp_list_path <- "/ccs/home/arodriguez/med112/task0101113/batch/pre-gwas/scripts_for_run/SNP_list.txt"
EUR_exclude_path <- "/gpfs/alpine/med112/proj-shared/data/genetic_annotations/R4.exclude_Imp3MAC20.HARE.EUR.txt.gz"
EUR_include_path <- "/ccs/home/arodriguez/med112/task0101113/batch/pre-gwas/scripts_for_run/R4.include_Imp3MAC20.HARE.EUR.txt.gz"
#for i in  /gpfs/alpine/med112/proj-shared/task0101113/output/pheCodes/inputs/bgen/*.pvar; do cut -f 3 $i|grep -v "#" | grep -v ID;done > $snp_list_path

.libPaths( c( .libPaths(), "/ccs/home/arodriguez/med112/task0101113/YoungDae_work/latest_2022030318/R_libs/") )
library(data.table)
snp_list <- fread(snp_list_path, header=FALSE)
exclude_list <- fread(EUR_exclude_path, header=FALSE)
include_list <- snp_list[! snp_list$V1 %in% exclude_list$V1,]
fwrite(include_list, EUR_include_path,  quote=FALSE, row.names=FALSE, col.names=FALSE)
