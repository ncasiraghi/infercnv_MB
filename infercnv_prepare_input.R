library( Matrix )

patients_folder <- list.files("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10XRNA5P",full.names = T)

patients <- c()
for(p in patients_folder){
  message(p)
  patients <- c(patients, grep(list.files(p,full.names = T),pattern = '-Nuclei|-PDX',value = TRUE))
}

main_outdir <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/rna/infercnv"

wd_list <- file.path(main_outdir,basename(patients))

for(wd in wd_list){
  
  dir.create(wd)
  
  message(wd)
  setwd(wd)
  
  cellranger_outs_folder_negative <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10XRNA5P/STP/STP-Nuclei/outs/"
  cell.annotation.name <- paste0("cellAnnotations_",basename(wd),".txt")
  counts.sparse.matrix.name <- paste0("sc.10x.counts_",basename(wd),".RData")
  
  # get count matrices from which tumor barcodes will be extracted
  
  cellranger_all_results <- grep(patients,pattern = basename(wd),value = T)
  
  cellranger_outs_folder_positive <- file.path(cellranger_all_results,"outs")
  
  # Accomodating 10X Data
  
  getCountMat <- function(cellranger_outs_folder){
    matrix_dir = file.path(cellranger_outs_folder,"filtered_feature_bc_matrix")
    barcode.path <- file.path(matrix_dir, "barcodes.tsv.gz")
    features.path <- file.path(matrix_dir, "features.tsv.gz")
    matrix.path <- file.path(matrix_dir, "matrix.mtx.gz")
    mat <- readMM(file = matrix.path)
    feature.names = read.delim(features.path,
                               header = FALSE,
                               stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path, 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    colnames(mat) = barcode.names$V1
    rownames(mat) = feature.names$V1
    return( mat )
  }
  
  getBarcodeNames <- function(cellranger_outs_folder){
    matrix_dir = file.path(cellranger_outs_folder,"filtered_feature_bc_matrix")
    barcode.path <- file.path(matrix_dir, "barcodes.tsv.gz")
    barcode.names = read.delim(barcode.path, 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    return( as.character(barcode.names[,1]) )
  }
  
  # Get count matrix
  
  mat.pos.list <- lapply( cellranger_outs_folder_positive, getCountMat )

  mat.neg.list <- lapply(cellranger_outs_folder_negative, getCountMat)
  
  n <- read.delim('/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/rna/infercnv/normalcells_barcodes/BarcodesALL.clean.tsv',stringsAsFactors = F)
  n <- n[which(n[,2] == 'sample1'),]
  m <- mat.neg.list[[1]]
  m <- m[,unique(n[,1])]
  mat.neg.list[[1]] <- m
  
  # add label for cells from the same sample
  for(i in seq(length(mat.pos.list))){
    colnames(mat.pos.list[[i]]) <- paste0(colnames(mat.pos.list[[i]]), "pos")
  }
  
  for(i in seq(length(mat.neg.list))){
    colnames(mat.neg.list[[i]]) <- paste0(colnames(mat.neg.list[[i]]), "neg")
  }
  
  check_and_merge <- function(sparse_matrix_list){
    if(length(sparse_matrix_list) < 2){
      return(sparse_matrix_list[[1]])
    }
    combtab <- combn(x = seq(length(sparse_matrix_list)),m = 2)
    
    genes_in_common <- list()
    for(i in seq(ncol(combtab))){
      genes_in_common[[i]] <- intersect(rownames(sparse_matrix_list[[combtab[1,i]]]), rownames(sparse_matrix_list[[combtab[2,i]]]))
    }
    
    cg <- genes_in_common[[1]]
    if(length(genes_in_common) > 1 ){
      for(i in 2:length(genes_in_common)){
        cg <- intersect(cg,genes_in_common[[i]])
      }
    }
    
    filter_genes <- function(x, genes){
      x <- x[genes,]
      return(x)
    }
    
    sparse_matrix_list <- lapply(sparse_matrix_list, filter_genes, cg)
    
    check.cols <- c()
    for(i in seq(ncol(combtab))){
      check.cols <- c(check.cols, identical(rownames(sparse_matrix_list[[combtab[1,i]]]), rownames(sparse_matrix_list[[combtab[2,i]]])) )
    }
    if(all(check.cols)){
      outmat <- do.call(cbind, sparse_matrix_list)
      return(outmat)
    } else {
      return(NULL)
    }
  }
  
  mat.pos <- check_and_merge(sparse_matrix_list = mat.pos.list)
  mat.neg <- check_and_merge(sparse_matrix_list = mat.neg.list)
  
  mat.matrix.sparse <- check_and_merge(sparse_matrix_list = list(mat.pos,mat.neg))
  
  dim(mat.pos)
  dim(mat.neg)
  dim(mat.matrix.sparse)
  
  # merge pos and neg
  save(mat.matrix.sparse, file = counts.sparse.matrix.name,compress = T)
  
  # cell annotations
  
  cellAnnotations <- data.frame(cell_ids=colnames(mat.matrix.sparse),
                                annotation=c(rep("tumor",as.numeric(lapply(mat.pos.list,ncol))),
                                             rep("normal",as.numeric(lapply(mat.neg.list,ncol)))
                                ),
                                stringsAsFactors = F)
  
  write.table(cellAnnotations, file=cell.annotation.name, col.names = F, row.names = F,quote=F, sep="\t")                   
  
  # genes ordering
  if( TRUE ){
    gene_ordering_file <- read.delim(file = "/icgc/dkfzlsdf/analysis/B260/users/n790i/tools/binning_the_genome/humangenes_biomart_GRCh37p13.sort.bed",stringsAsFactors = F,skip = 1,header = F)
    gene_ordering_file <- gene_ordering_file[which(gene_ordering_file$V1 != "MT"),]
    gene_ordering_file$index <- gene_ordering_file$V1
    
    gene_ordering_file$index[which(gene_ordering_file$index == "X")] <- 23
    gene_ordering_file$index[which(gene_ordering_file$index == "Y")] <- 24
    gene_ordering_file$index <- as.numeric(gene_ordering_file$index)
    gene_ordering_file$V1 <- paste0("chr",gene_ordering_file$V1)
    gene_ordering_file <- gene_ordering_file[with(gene_ordering_file,order(index)),]
    write.table(gene_ordering_file[,c(5,1:3)], file='gene_ordering_file.txt', col.names = F, row.names = F,quote=F, sep="\t")
  }
  
}
