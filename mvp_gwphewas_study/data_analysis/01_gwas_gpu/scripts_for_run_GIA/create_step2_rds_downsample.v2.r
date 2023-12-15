#!/usr/bin/env Rscript
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

# script will create the df for a particular phenotype and group

colClasses = c("character", "character", "character", "character", "character", "character", "character", "character", "character", "character", "character", "character", "character", "character")
col.names = c("--bgenFile=", "--bgenFileIndex=", "--sampleFile=", "--GMMATmodelFile=", "--varianceRatioFile=", "--SAIGEOutputFile=", "--minMAF=", "--minMAC=", "--numLinesOutput=", "--IsOutputNinCaseCtrl=", "--IsOutputAFinCaseCtrl=", "--LOCO=", "--idstoIncludeFile", ">> ")
df <- read.table(text = "", colClasses = colClasses, col.names = col.names)
df_orig <- read.table(text = "", colClasses = colClasses, col.names = col.names)

phenoID = args[1]
group = args[2]
outrds = args[3]
chrms = unlist(strsplit(args[4],",",fixed=T))
#chrms = paste("chr", chrms, "\\.", sep="")
#step1_output = args[4]
include_file = args[5]

sampF="Release4.dose.mvpcoreid.psam"
rdaF= paste("downsample.phewas.ld.maf", phenoID, group, "out.rda", sep=".")
varF= paste("downsample.phewas.ld.maf", phenoID, group, "out.varianceRatio.txt", sep=".")

orig_rdaF= paste("phewas.ld.maf", phenoID, group, "out.rda", sep=".")
orig_varF= paste("phewas.ld.maf", phenoID, group, "out.varianceRatio.txt", sep=".")

base_path = "/gpfs/alpine/med112/proj-shared/task0101113/output/GIA_ANC_Run"
step1_input_path= paste(base_path, phenoID, group, "step1", sep="/")
output_path= paste(base_path, phenoID, group, "step2", sep="/")
bgen_path="/gpfs/alpine/med112/proj-shared/task0101113/output/pheCodes/inputs/bgen"
bgen_col = paste("", list.files(bgen_path, pattern="bgen$", full.names=TRUE), sep="")
bgen_index_col = paste("", paste(bgen_col, "bgi", sep="."), sep="")
samp_col = paste("", replicate(length(bgen_col), paste(bgen_path, sampF, sep="/")), sep="")
GMMAT_col = paste("", replicate(length(bgen_col),paste(step1_input_path, rdaF, sep="/")), sep="")
var_col = paste("", replicate(length(bgen_col),paste(step1_input_path, varF, sep="/")), sep="")
orig_GMMAT_col = paste("", replicate(length(bgen_col),paste(step1_input_path, orig_rdaF, sep="/")), sep="")
orig_var_col = paste("", replicate(length(bgen_col),paste(step1_input_path, orig_varF, sep="/")), sep="")
output_col = paste("", paste(output_path, paste(list.files(bgen_path, pattern="bgen$"), phenoID, group, "downsample.assoc.txt", sep="."), sep="/"), sep="")
orig_output_col = paste("", paste(output_path, paste(list.files(bgen_path, pattern="bgen$"), phenoID, group, "downsample_orig.assoc.txt", sep="."), sep="/"), sep="")
minMAF_col =  paste("", replicate(length(bgen_col), 0.0001), sep="")
minMAC_col = paste("", replicate(length(bgen_col), 1), sep="")
numLinesOutput_col = paste("", replicate(length(bgen_col), 10000), sep="")
IsOutputNinCaseCtrl_col = paste("", replicate(length(bgen_col), "TRUE"), sep="")
IsOutputAFinCaseCtrl_col = paste("", replicate(length(bgen_col), "TRUE"), sep="")
LOCO_col = paste("", replicate(length(bgen_col), FALSE), sep="")
include_col= paste("", replicate(length(bgen_col), include_file), sep="")
outputlog_col = paste("", paste(output_path, paste(list.files(bgen_path, pattern="bgen$"), phenoID, group, "downsample.log", sep="."), sep="/"), sep=" ")
orig_outputlog_col = paste("", paste(output_path, paste(list.files(bgen_path, pattern="bgen$"), phenoID, group, "downsample_orig.log", sep="."), sep="/"), sep=" ")

df <- data.frame(bgen_col, bgen_index_col, samp_col, GMMAT_col, var_col, output_col, minMAF_col, minMAC_col, numLinesOutput_col, IsOutputNinCaseCtrl_col, IsOutputAFinCaseCtrl_col, LOCO_col, outputlog_col, include_col)
df <- df[grepl(paste(chrms, collapse= "|"), df$outputlog_col, fixed = FALSE) ,]
df <- as.matrix(df)
saveRDS(object=df, file=outrds)

df_orig <- data.frame(bgen_col, bgen_index_col, samp_col, orig_GMMAT_col, orig_var_col, orig_output_col, minMAF_col, minMAC_col, numLinesOutput_col, IsOutputNinCaseCtrl_col, IsOutputAFinCaseCtrl_col, LOCO_col, orig_outputlog_col, include_col)
df_orig <- df_orig[grepl(paste(chrms, collapse= "|"), df_orig$orig_outputlog_col, fixed = FALSE) ,]
df_orig <- as.matrix(df_orig)
outrds_orig = paste(outrds, "original.rds", sep=".")
saveRDS(object=df_orig, file=outrds_orig)
