#!/bin/sh 
#BSUB -J infercnv
#BSUB -q verylong 
#BSUB -e bsub_infercnv.log 
#BSUB -o bsub_infercnv.txt 

module load R/3.6.0
#Rscript /icgc/dkfzlsdf/analysis/B260/users/n790i/infercnv_MB/infercnv_analysis.R /icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/rna/infercnv/ST1R-Nuclei
Rscript /icgc/dkfzlsdf/analysis/B260/users/n790i/infercnv_MB/infercnv_analysis.R /icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/rna/infercnv/4M67-Nuclei
Rscript /icgc/dkfzlsdf/analysis/B260/users/n790i/infercnv_MB/infercnv_analysis.R /icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/rna/infercnv/BT084-PDX
Rscript /icgc/dkfzlsdf/analysis/B260/users/n790i/infercnv_MB/infercnv_analysis.R /icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/rna/infercnv/MB243-Nuclei
Rscript /icgc/dkfzlsdf/analysis/B260/users/n790i/infercnv_MB/infercnv_analysis.R /icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/rna/infercnv/RCMB18-PDX
Rscript /icgc/dkfzlsdf/analysis/B260/users/n790i/infercnv_MB/infercnv_analysis.R /icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/rna/infercnv/ST1R-PDX
Rscript /icgc/dkfzlsdf/analysis/B260/users/n790i/infercnv_MB/infercnv_analysis.R /icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/rna/infercnv/STP-Nuclei
Rscript /icgc/dkfzlsdf/analysis/B260/users/n790i/infercnv_MB/infercnv_analysis.R /icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/rna/infercnv/STP-PDX







