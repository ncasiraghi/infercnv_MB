#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=1){
  message("\nError!\nusage: Rscript infercnv_output_matrix.R /path/to/infercnv/results\n")
  quit()
}

wd <- args[1]
input <- file.path(wd,'results','infercnv.observations.txt')
finalmat <- read.table(file = input,header = T,row.names = 1)
colnames(finalmat) <- gsub(colnames(finalmat),pattern = '\\.1pos',replacement = '')
save(finalmat,file = file.path(wd,paste0(basename(wd),'_infercnv.matrix.RData')),compress = TRUE)
