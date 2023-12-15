# EUR:
cd .//output/GIA_ANC_Run/SUBMIT
for i in ../*/EUR; do echo $i; find $i/step2/*.txt -type f -print0 | xargs -0 -L1 bash -c 'test "$(tail -c 1 "$0")" && echo "No new line at end of $0"'; done > incomplete_EUR_chunks.txt
for i in ../*/EUR/step2; do echo $i; grep -L "Analysis done" $i/*.log; done > incomplete_eur_chunks.txt
for i in `grep log incomplete_eur_chunks.txt`; do echo $i; dir=`dirname $i`; base=`basename $i .log`; txt=$dir/$base; echo $txt; rm $txt; rm $i; done
for i in `grep  "No new line" incomplete_EUR_chunks.txt |  grep EUR | cut -f 7 -d " " | more`; do echo $i; rm $i; rm $i.log;done

Rscript ../../../batch/pre-gwas/gwPheWAS-Summit/scripts_for_run_GIA/merge_step2_rds_gpu_v1.r EUR
for i in ../*/EUR/step2; do size=`ls $i/*.gpu.txt.log| wc -l`; echo $i $size;done > EUR.counts.completed.txt


# in R
x_eur <- readRDS(".//output/GIA_ANC_Run/SUBMIT/complete_step2.finishedstep1tasks.gpu.EUR.rds")
eur_phecodes <- read.csv(".//output/GIA_ANC_Run/SUBMIT/EUR.counts.completed.txt", sep=" ", header=F)
# play around with the numbers to figure out which ones need to be done
incomplete_eur_phecodes <- do.call(rbind, strsplit(eur_phecodes[eur_phecodes$V2 < 219 & eur_phecodes$V2 > 0,]$V1, '/'))[,2]
x_eur_left <- x_eur[grepl("Platelet_Mean_INT/ASN",x_eur$outputlog_col),]
for (i in incomplete_eur_phecodes) { x_eur_left <- rbind(x_eur_left, x_eur[grepl(paste(i,"EUR", sep="/"),x_eur$outputlog_col),])}
x_eur_left <- unique(x_eur_left)
saveRDS(x_eur_left, ".//output/GIA_ANC_Run/SUBMIT/complete_step2.leftover.finishedstep1tasks.gpu.EUR.rds")

# then submit the leftovers
