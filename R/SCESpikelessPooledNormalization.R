
#### Load Libraries
library('scran')
library('FactoMineR')
library('Rtsne')

##### Data pre-processing ==================================================================================================
##### Reading and ordering non auxillary tables------------------------------------------

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



#### Plotting functions =====================================================================================================================
  ## Plot PCA
plot.pca<-function(pca, title, color="black"){
  plot(pca$ind$coord[,1:2], 
       pch=19,
       xlab=paste0("PC 1\n(", round(pca$eig[1,2], digits=2), "%)"),
       ylab=paste0("PC 2 (", round(pca$eig[2,2], digits=2), "%)"),
       main=title,
       col=color)
  legend("topright",
         legend=as.character(levels(color)), 
         pch=rep(19, length(levels(color))), 
         col=as.factor(levels(color)), 
         bty="n",
         xpd=NA)
}  
  ## Plot TSNE
plot.tsne<-function(tsne, title, color="black"){
  plot(tsne$Y, 
       pch=19,
       main=title,
       col=color,
       xlab="Dim1",
       ylab="Dim2")
  legend("topright",
         legend=as.character(levels(color)), 
         pch=rep(19, length(levels(color))), 
         col=as.factor(levels(color)), 
         bty="n",
         xpd=NA)
}

  ## Realise all plots for data
plot.data<-function(matrix, color="black", normalized=FALSE) {
  
  # Set Data name
  if(normalized) data<-"Normalized Data"
  else data<-"Raw Data"
  
  # Plot unscaled PCA
  pca<-PCA(t(log(matrix+1)), graph=FALSE, scale.unit=FALSE)
  pdf(paste0("Scran_", sub(" ", "", data), "_unscaledPCA.pdf"))
  plot.pca(pca, color=color, title=paste0(data, ' PCA (unscaled, log transformed)'))
  dev.off()
  
  # Plot scaled PCA
  pca<-PCA(t(matrix), graph=FALSE, scale.unit=TRUE)
  pdf(paste0("Scran_", sub(" ", "", data), "_scaledPCA.pdf"))
  plot.pca(pca, color=color, title=paste0(data, ' PCA (scaled, log transformed)'))
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
  pdf(file=paste0("Scran_", sub(" ", "", data), "_tsne.pdf"), width=8*length(n.neighbors), height=8)
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
               color.by="Condition", preclustering=FALSE, 
               min_cluster_size=100, force.positive=FALSE,
               detection = 1, nCells = 1,
               cell_cycle=FALSE, organism='mus_musculus',
               output){
  
  ### This function integrates all the parameters given by the user and write the normalized 
  ### expression matrix and the associated size factor in the cells phenotype file.
  
  # Importing Data
  sce <- import_SCE(file)
  
  
  if(cell_cycle){
    if(organism=='mus_musculus'){
      pairs<-readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
    }else{
      if(organism =='homo_sapiens'){
        pairs<-readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
      }else{
        stop ("Organism for cell cycle annotation isn't valid.")
      }
    }
    assigned<-cyclone(sce, pairs) 
    SummarizedExperiment::colData(sce)$Phase <- assigned$phases
  }
  
  
  # Treating expression matrix
  fsce <- filter_SCE(sce, mode)
  
    # Filter genes on expression
  kept<- rowSums(count >= detection) >= nCells
  if (length(kept) == 0) {
    warning("No genes passed expression criterion, continuing with whole dataset")
  } else {
    fsce <- fsce[kept,]
  } 
   
  # Plotting raw data
  plot.data( SummarizedExperiment::assay(sce, "counts"), 
            color=as.factor(pheno[,grep(color.by, colnames(pheno))]), 
            normalized=FALSE)
  
 
  
  # Calculate pooling sizes vector
  pooling_size s<- floor(ncol(sce)/20)*10
  pooling_sizes <- seq(from=pooling_sizes/5, to=pooling_sizes, by=pooling_sizes/5)
  
  # Calculating size factors
    if(preclustering){
      clusters <- quickCluster(fsce, min.size=min_cluster_size)
      fsce <- computeSumFactors(fsce, cluster=clusters, sizes=pooling_sizes/length(levels(clusters)), positive=force.positive, get.spikes=FALSE, errors=FALSE)
    }else{
      fsce <- computeSumFactors(fsce, cluster=NULL, sizes=pooling_sizes, positive=force.positive, get.spikes=FALSE, errors=FALSE)
    }
  BiocGenerics::sizeFactors(sce) <- BiocGenerics::sizeFactors(fsce)
  rm(fsce)
  
  # Normalizing data
  sce <- BiocGenerics::normalize(sce, return_log = FALSE, log_exprs_offset = 0)
  
    # Annotating for cell cycle
 
  saveRDS(sce, output1)
  
  
  # Plotting Normalized Data
  plot.data(SummarizedExperiment::assay(sce, "normcounts"), 
            color=as.factor(pheno[,grep(color.by, colnames(pheno))]), 
            normalized=TRUE)
  
  # Clean directory
  if(file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
  
  # End
  return("Normalization process success")
}
  
        ### Launch commands
args<-commandArgs(TRUE)
main(args[1], args[2], args[3],
      as.logical(args[4]), as.numeric(args[5]), as.logical(args[6]),
      as.numeric(args[7]), as.numeric(args[8]),
      as.logical(args[9]), args[10],
      args[11])
