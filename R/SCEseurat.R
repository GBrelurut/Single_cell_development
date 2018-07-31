library(SingleCellExperiment)

##### Data pre-processing ===============================================================

#### Testing parameters

as.valid.logical<-function(test) {
  res<-as.logical(test)
  if(is.na(res)) stop("Logical parameters must be : one of true, TRUE, True, False, false, or FALSE")
  return(res)
}


##### Reading and ordering non auxillary tables------------------------------------------

## Import SCE
import_SCE <- function(file) {
  return(readRDS(file))
}


#### Apply corrections =================================================

normalizeSeurat<-function(object) {
 object@data<-log(t(t(object@raw.data)/object@meta.data$SizeFactors) + 1)
  return(object)
}

#### HVG detection =====================================================

### Return highly variable genes identify by scran
hvg.with.scran<-function(seurat.object, spike = TRUE, 
                         low.mean = 0.01, high.mean = 5, bio = 1 ) {
  
  # construct object
  SCE <- Seurat::as.SingleCellExperimen(seurat.object)
  
  # find spikes if needed (can be long)
  if(spike && any(grepl("spike", rowData(SCE)$Type))) {
    isSpike(SCE, "ERCC") <- rowData(SCE)$Type == "spike"
  }
  
  # find HVGs
  fit <- scran::trendVar(SCE, trend = "loess",  use.spikes = (! is.null(isSpike(SCE))) )
  dc <- scran::decomposeVar(SCE, fit)
  
  # Plot results
  pdf("scranFittedVariability.pdf")
  plot(dc$mean, dc$total, xlab="Mean log-expression", 
         ylab="Variance", pch ="19", col = "gray")
  o <- order(dc$mean)
  lines(dc$mean[o], dc$tech[o], col="blue", lwd=2)
  if (! is.null(isSpike(SCE))) points(fit$mean, fit$var, col="blue", pch=19)
  dev.off()
  
  # get and return hvg
  hvg <- rownames(dc)[which( dc$mean > low.mean && dc$mean <= high.mean && dc$bio >= bio)]
  return(hvg)
}


#### Cluster quality metrics ===========================================

#### Silhouette plot ---------------------------------------------------

### Calculate and plot silhouette from seurat object with clusters
### parameters :
###   seurat.object : seurat object to use
###   pc.used : principal component used for clustering
###   distance : distance to use for silhouette calculation
###   plot : weither to plot silhouette or not
### returns : silhouette object from cluster package

silhouette.seurat<-function(seurat.object, pc.used, distance="euclidean", plot=TRUE) {
  library(cluster)
  coord <- seurat.object@dr$pca@cell.embeddings
  coord <- coord[,pc.used]
  clusters <- seurat.object@ident
  d <- dist(coord, method=distance)
  sil<-silhouette(as.numeric(clusters), dist=d)
  if(plot){
    par(mfrow = c(1, 1))
    plot(sil, col=as.factor(clusters[order(clusters, decreasing=FALSE)]), main="Silhouette plot of Seurat clustering")
    abline(v=mean(sil[,3]), col="red4", lty=2)
  }
  return(sil)
}


#### Main function 
main<-function(file,  # input file
               normalize = TRUE, scale = TRUE, # normalization and scaling
               n.cells = 0, detection = 10, # filtering genes on detection
               hvg.detection = "none", low.mean, high.mean, var.threshold, spike = FALSE, # filtering on variability
               n.cores = 2, jackstraw.replicate = 1000, jackstraw.prop = 0.1, sig.threshold = 0.05, score.threshold = 1e-5, # selecting significative dimensions
               resolution = 0.6, k.param = 30, algorithm="Louvain", sparse = FALSE # clustering options
               ) {
  
  library(methods)
  
  # Checking parameters
  if( ! hvg.detection %in% c("scran", "seurat", "none")) stop("hvg.detection must be one of : scran, seurat, or none")
  if( ! algorithm %in% c("Louvain", "Louvain.multilevel", "SLM")) stop("algorithm must be one of : Louvain, Louvain.multilevel, or SLM")
  
  ## import data and construct Seurat object
  object <- import_SCE(file)
  SF <- BiocGenerics::sizeFactors(object)
  names(SF) <- colnames(object)
  assay(object, "logcounts") <- log(counts(object) + 1) 
  object <- Seurat::as.seurat(object)
  object <- Seurat::AddMetaData(object, SF, "SizeFactors")
  rm(SF)
  
  ## Normalize data
  if(normalize) object <- normalize.seurat(object)
  else object@data <- log( object@raw.data + 1 )
 
  ## Filter genes on expression
  if(n.cells) {
    kept <- which(rowSums( object@raw.data >= detection) >=  n.cells )
    object@data <- object@data[kept,]
  } 
  
  ## Scale data
  if(scale) object@scale.data <- t(scale(t(object@data)))
  else object@scale.data <- object@data
 
  ## Select HVG
  genes<-import_controls(genes.file)
  if( hvg.detection == "none") hvg <- rownames(object@scale.data)
  if( hvg.detection == "scran") hvg <- hvg.with.scran(object, genes, spike, low.mean, 
                                                      high.mean, bio = var.threshold ) 
  if( hvg.detection == "seurat") {
    pdf("seuratFittedVariation.pdf")            # Default parameters
    object <- Seurat::FindVariableGenes(object, # mean.function = ExpMean, dispersion.function = LogVMR,
                          do.plot = TRUE,
                          x.low.cutoff = low.mean, x.high.cutoff = high.mean,
                          y.cutoff = var.threshold)
    dev.off()
    hvg <- object@var.genes
  } 
    ## If no genes are kept keep all genes
  if(length(hvg) == 0) {
    hvg <- rownames(object@scale.data)
    warning("No Highly variable gene found. Continue with all genes.")
  }
  
  ## Save genes used for dimensionnality reduction
  object@var.genes <- hvg 
  
  ## Perform PCA and JackStraw
  object <- Seurat::RunPCA(object, pc.genes = hvg, do.print = FALSE, 
                pcs.compute = min(ncol(object@scale.data), length(hvg))/2 - 1) # Optimization prevents to compute 50% of eigen values or more
  parallelized <- FALSE
  if(n.cores > 1) parallelized <- TRUE
  object <- Seurat::JackStraw(object, num.pc = ncol(object@dr$pca@gene.loadings), prop.freq = jackstraw.prop,
                      num.replicate = jackstraw.replicate, do.par = parallelized, num.cores = n.cores)
  
  ## Select significative PCs
  
  pc.scores <- object@dr$pca@jackstraw@emperical.p.value
  p.value <- 0
  i <- 1
  pcs <-NULL
  while(  i < ncol(pc.scores)) {
    n.obs <- sum(pc.scores[,i] <= score.threshold)
    if(n.obs == 0 ) {
      p.value <- 1
    }
    tot <- nrow(pc.scores)
    n.theo <- tot * score.threshold
    p.value <- prop.test(c(n.obs, n.theo), c(tot, tot))$p.val
    if(p.value <= sig.threshold) pcs <- c(pcs, i)
    i <- i + 1
  }
  rm(pc.scores)
  
  # Treat few pcs 
  if(length(pcs) < 2) pcs <- 1:2
  
  ## plot and save results
  object <- Seurat::JackStrawPlot(object, PCs = 1:max(pcs + 1), score.thresh = score.threshold)
  ggplot2::ggsave("PCSignificativityPlots.pdf", object@dr$pca@misc$jackstraw.plot , device = "pdf")
  
  ## Produce PCs heatmaps
  ncells<-min(200, ncol(object@scale.data))
  pdf("PCsHeatmap.pdf")
  Seurat::PCHeatmap(object, pc.use = pcs, cells.use = ncells, 
            do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
  dev.off()
  
  ## Cluster cells
  if(algorithm == "Louvain") algorithm <- 1
  if(algorithm == "Louvain.multilevel") algorithm <- 2
  if(algorithm == "SLM") algorithm <- 3
  
  object <- Seurat::FindClusters(object, dims.use = pcs, resolution = resolution,
                         print.output = TRUE, save.SNN = TRUE,  k.param = k.param, 
                         algorithm = algorithm, do.sparse = sparse, temp.file.location=getwd())
  
  ## Visualize using TSNE
  object <- Seurat::RunTSNE(object, dims.use = pcs, do.fast=FALSE)
  p <- Seurat::TSNEPlot(object, do.return = TRUE)
  ggplot2::ggsave("clusterPlot.pdf", plot = p, device = "pdf")
  
  ## Evaluate with silhouette
  pdf("silhouettePlot.pdf")
  sil <- silhouette.seurat(object, pc.used = pcs, distance="euclidean", plot=TRUE)
  dev.off()
  
  ## Save clustering
  saveRDS(object, file = "seurat.RDS")
}

## Treat Parameters
args <- commandArgs(TRUE)
arg2 <- as.valid.logical(args[2])
arg3 <- as.valid.logical(args[3])
arg10 <- as.valid.logical(args[10])
arg19 <- as.valid.logical(args[19])

## Launch main function
main(args[1], # input file
     arg2, arg3, # correction options
     as.numeric(args[4]), as.numeric(args[5]), # filtering genes on detection
     args[6], as.numeric(args[7]), as.numeric(args[8]), as.numeric(args[9]), arg10, #filtering on variability
     as.numeric(args[11]), as.numeric(args[12]), as.numeric(args[13]), as.numeric(args[14]), as.numeric(args[15]), # jackstraw options
     as.numeric(args[16]), as.numeric(args[17]), args[18], arg19 # clustering options
     )
