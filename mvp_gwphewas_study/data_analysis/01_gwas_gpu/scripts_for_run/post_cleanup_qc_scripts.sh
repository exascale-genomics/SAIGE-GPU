# AFR:
cd .//output/HARE_ANC_Run/SUBMIT
for i in ../*/AFR; do echo $i; find $i/step2/*.txt -type f -print0 | xargs -0 -L1 bash -c 'test "$(tail -c 1 "$0")" && echo "No new line at end of $0"'; done > incomplete_AFR_chunks.txt
for i in ../*/AFR/step2; do echo $i; grep -L "Analysis done" $i/*.log; done > incomplete_afr_chunks.txt
for i in `grep log incomplete_afr_chunks.txt`; do echo $i; dir=`dirname $i`; base=`basename $i .log`; txt=$dir/$base; echo $txt; rm $txt; rm $i; done
for i in `grep  "No new line" incomplete_AFR_chunks.txt |  grep AFR | cut -f 7 -d " " | more`; do echo $i; rm $i; rm $i.log;done

Rscript ../../../batch/pre-gwas/gwPheWAS-Summit/scripts_for_run/merge_step2_rds_gpu_v1.r AFR
for i in ../*/AFR/step2; do size=`ls $i/*.gpu.txt.log| wc -l`; echo $i $size;done > AFR.counts.completed.txt


# in R
x_afr <- readRDS(".//output/HARE_ANC_Run/SUBMIT/complete_step2.finishedstep1tasks.gpu.AFR.rds")
afr_phecodes <- read.csv(".//output/HARE_ANC_Run/SUBMIT/AFR.counts.completed.txt", sep=" ", header=F)
# play around with the numbers to figure out which ones need to be done
incomplete_afr_phecodes <- do.call(rbind, strsplit(afr_phecodes[afr_phecodes$V2 < 219 & afr_phecodes$V2 > 110,]$V1, '/'))[,2]
x_afr_left <- x_afr[grepl("Platelet_Mean_INT/ASN",x_afr$outputlog_col),]
for (i in incomplete_afr_phecodes) { x_afr_left <- rbind(x_afr_left, x_afr[grepl(paste(i,"AFR", sep="/"),x_afr$outputlog_col),])}
x_afr_left <- unique(x_afr_left)
saveRDS(x_afr_left, ".//output/HARE_ANC_Run/SUBMIT/complete_step2.leftover.finishedstep1tasks.gpu.AFR.rds")

# then submit the leftovers
