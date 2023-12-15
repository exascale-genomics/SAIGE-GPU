#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# script will create the df for a particular phenotype and group

colClasses = c("character", "character", "character", "character", "character", "character", "character", "character", "character", "character", "character", "character", "character")
col.names = c("--bgenFile=", "--bgenFileIndex=", "--sampleFile=", "--GMMATmodelFile=", "--varianceRatioFile=", "--SAIGEOutputFile=", "--minMAF=", "--minMAC=", "--numLinesOutput=", "--IsOutputNinCaseCtrl=", "--IsOutputAFinCaseCtrl=", "--LOCO=", ">> ")
df <- read.table(text = "", colClasses = colClasses, col.names = col.names)

phenoID = args[1]
group = args[2]
outrds = args[3]
#step1_output = args[4]

sampF="Release4.dose.mvpcoreid.psam"
rdaF= paste("phewas.ld.maf", phenoID, group, "out.rda", sep=".")
varF= paste("phewas.ld.maf", phenoID, group, "out.varianceRatio.txt", sep=".")
base_path = "./output/GIA_ANC_Run"
step1_input_path= paste(base_path, phenoID, group, "step1", sep="/")
output_path= paste(base_path, phenoID, group, "step2", sep="/")
bgen_path="./output/pheCodes/inputs/bgen"

bgen_col = paste("", list.files(bgen_path, pattern="bgen$", full.names=TRUE), sep="")
bgen_index_col = paste("", paste(bgen_col, "bgi", sep="."), sep="")
samp_col = paste("", replicate(length(bgen_col), paste(bgen_path, sampF, sep="/")), sep="")
GMMAT_col = paste("", replicate(length(bgen_col),paste(step1_input_path, rdaF, sep="/")), sep="")
var_col = paste("", replicate(length(bgen_col),paste(step1_input_path, varF, sep="/")), sep="")
output_col = paste("", paste(output_path, paste(list.files(bgen_path, pattern="bgen$"), phenoID, group, "assoc.txt", sep="."), sep="/"), sep="")
minMAF_col =  paste("", replicate(length(bgen_col), 0.0001), sep="")
minMAC_col = paste("", replicate(length(bgen_col), 1), sep="")
numLinesOutput_col = paste("", replicate(length(bgen_col), 10000), sep="")
IsOutputNinCaseCtrl_col = paste("", replicate(length(bgen_col), "TRUE"), sep="")
IsOutputAFinCaseCtrl_col = paste("", replicate(length(bgen_col), "TRUE"), sep="")
LOCO_col = paste("", replicate(length(bgen_col), FALSE), sep="")
outputlog_col = paste("", paste(output_path, paste(list.files(bgen_path, pattern="bgen$"), phenoID, group, "log", sep="."), sep="/"), sep=" ")

df <- rbind(df, as.matrix(data.frame(bgen_col, bgen_index_col, samp_col, GMMAT_col, var_col, output_col, minMAF_col, minMAC_col, numLinesOutput_col, IsOutputNinCaseCtrl_col, IsOutputAFinCaseCtrl_col, LOCO_col, outputlog_col)))
df <- as.matrix(df)
saveRDS(object=df, file=outrds)

