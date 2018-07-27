#! /usr/bin/env Rscript

#### Load Libraries
library('FactoMineR')
library('Rtsne')
library('SingleCellExperiment')

##### Data pre-processing ==================================================================================================
  ##### Reading and ordering non auxillary tables

## Import Data
import_SCE <- function(file) {
  return(readRDS(file))
}

#### Processing Expression Matrix -----------------------------------------------------

filter_SCE <- function(SCE, mode = "Nuclear") {
  
  # Check mode	
  if(! mode %in% c("Nuclear", "Endogenous", "All")) stop("Matrix filtering should be one of All, Nuclear or Endogenous")
  
  # If All filter nothing
  if(mode == "All") return(SCE)
  
  # Get and remove spikes
  spikes <- SummarizedExperiment::rowData(SCE)$Type == "spike"
  SCE <- SCE[! spike,]
  
  if(mode == "Endogenous") return(SCE)
  
  # Get and remove mitochondrial genes
  mt <- SummarizedExperiment::rowData(SCE)$Type == "mitochondrial"
  return(SCE[! mt,])
  
}



## Function for length correction
correct_gene_length<-function(sce) {
  lengths<-rowData(sce)$Length
  sce@assays <- SummarizedExperiment::Assays(list(counts = SummarizedExperimen::assay(sce, "counts"), corr_counts = SummarizedExperimen::assay(sce, "counts")/lengths))
  return(sce)
}


## SF function
calculate_SF<-function(sce, mode='CPM') {
  if(mode == 'TPM' && ! SummarizedExperiment::rowData(sce)$GeneID ) {
    sce <- correct_gene_length(sce)
    BiocGenerics::sizeFactors(sce) <- colSums(SummmarizedExperiment::assay(sce, "corr_counts"))/1e6
  } else {
    BiocGenerics::sizeFactors(sce) <- colSums(SummmarizedExperiment::assay(sce, "counts"))/1e6
  }
  return(BiocGenerics::normalize(sce, return_log = FALSE, log_exprs_offset = 0))
}

#### Plotting functions =====================================================================================================================
  ## Plot PCA
plot.pca<-function(pca, title, color='black'){
  plot(pca$ind$coord[,1:2], 
       pch=19,
       xlab=paste0('PC 1\n(', round(pca$eig[1,2], digits=2), '%)'),
       ylab=paste0('PC 2 (', round(pca$eig[2,2], digits=2), '%)'),
       main=title,
       col=color)
  legend('topright',
         legend=as.character(levels(color)), 
         pch=rep(19, length(levels(color))), 
         col=as.factor(levels(color)), 
         bty='n',
         xpd=NA)
}  
  ## Plot TSNE
plot.tsne<-function(tsne, title, color='black'){
  plot(tsne$Y, 
       pch=19,
       main=title,
       col=color,
       xlab='Dim1',
       ylab='Dim2')
  legend('topright',
         legend=as.character(levels(color)), 
         pch=rep(19, length(levels(color))), 
         col=as.factor(levels(color)), 
         bty='n',
         xpd=NA)
}

  ## Realise all plots for data
plot.data<-function(matrix, color='black', normalized=FALSE) {
  
  # Set Data name
  if(normalized) data<-'Normalized Data'
  else data<-'Raw Data'
  
  # Plot unscaled PCA
  pca<-PCA(t(matrix), graph=FALSE, scale.unit=FALSE)
  pdf(paste0('Sum_', sub(' ', '', data), '_unscaledPCA.pdf'))
  plot.pca(pca, color=color, title=paste0(data, ' PCA (unscaled)'))
  dev.off()
  
  # Plot scaled PCA
  pca<-PCA(t(matrix), graph=FALSE, scale.unit=TRUE)
  pdf(paste0('Sum_', sub(' ', '', data), '_scaledPCA.pdf'))
  plot.pca(pca, color=color, title=paste0(data, ' PCA (scaled)'))
  dev.off()
  rm(pca)
  
  # Plot Tsne on log transformed data
    # Get perplexity max value
  n.neighbors<-floor(ncol(matrix)/3.5)
    # Calculate perplexity intermediate values
  n.neighbors<-c(2, min(floor(n.neighbors/5), 30), floor(n.neighbors/2), n.neighbors)
    # Clean vector
  n.neighbors<-unique(n.neighbors[which(n.neighbors>=2)])
    # Calculate and plot TSNE
  pdf(file=paste0('Sum_', sub(' ', '', data), '_tsne.pdf'), width=8*length(n.neighbors), height=8)
  layout(matrix(c(1:length(n.neighbors)), nrow=1, byrow=TRUE))
  for(nn in n.neighbors) {
    tsne<-Rtsne(t(log(matrix+1)), pca=FALSE, perplexity=nn)
    plot.tsne(tsne, color=color, title=paste0(data, ' T-SNE plot \n(log transformed, perplexity = ', nn, ')'))
  }
  dev.off()
  layout(matrix(c(1), nrow=1, byrow=TRUE))
}

##### Main function =============================================================================================================================================================================
main<-function(file, mode='Nuclear', 
               color.by='Condition', norm='CPM', 
               output){
  
  ### This function integrates all the parameters given by the user and write the normalized 
  ### expression matrix and the associated size factor in the cells phenotype file.
  
  # Importing Data
  sce <- import_SCE(file)
  
  # filtering genes
  fsce <- filter_SCE(sce, mode = mode )
  sce <- filter_SCE(sce, mode = "Endogenous")
  
  # Plotting raw data
  plot.data(SummarizedExperiment::assay(sce, "counts"), 
            color=as.factor(pheno[,grep(color.by, colnames(pheno))]), 
            normalized=FALSE)
  

  # Calculating size factors
  if(norm == 'TPM' && is.null(SummarizedExperiment::rowData(sce)$GeneID)) {
      warning('TPM normalization is not available for gene counting, shifting to CPM')
      norm <- 'CPM'
  }
  
  sce <- calculate_SF(sce, norm)
  
  # Saving
  saveRDS(sce, output)
 

  # Plotting Normalized Data
  plot.data(SummarizedExperiment::assay(sce, "normcounts"), 
            color=as.factor(pheno[,grep(color.by, colnames(pheno))]), 
            normalized=TRUE)
  
  # Clean directory
  if(file.exists('Rplots.pdf')) file.remove('Rplots.pdf')
  
}
  
        ### Launch commands
args<-commandArgs(TRUE)
main(args[1], args[2],
     args[3], args[4],
     args[5])
