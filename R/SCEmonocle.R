#! /usr/bin/env Rscript

## Import arguments
args<-commandArgs(TRUE)

## Load library
library('monocle')

### Functions ---------------------------------------------------------------------------
## Import Data
import_SCE <- function(file) {
  return(readRDS(file))
}


## Create Cell Data Set
cellDataSet<-function(sce, fam, detection.threshold) {
  ## Set option verctor
  dist.family<-c("negbinomial.size", "negbinomial", "tobit", "gaussianff") # Last families are untested, but conserved to help module update
  # Check parameter
  if(! fam %in% dist.family) stop(warning("Undifined distribution family"))
  # Create CDS
  if(fam == dist.family[1]) {
    CDS<-newCellDataSet(as.matrix(SummarizedExperiment::assay(sce, "counts")),
                        phenoData = data.frame(SummarizedExperiment::colData(sce)),
                        featureData = data.frame(SummarizedExperiment::rowData(sce)),
                        expressionFamily = negbinomial.size(), 
                        lowerDetectionLimit = detection.threshold)
  }
  
  if(fam == dist.family[2]) {
    CDS<-newCellDataSet(as.matrix(SummarizedExperiment::assay(sce, "counts")),
                        phenoData = data.frame(SummarizedExperiment::colData(sce)),
                        featureData = data.frame(SummarizedExperiment::rowData(sce)), 
                        expressionFamily = negbinomial(), 
                        lowerDetectionLimit = detection.threshold)
  }
  
  if(fam == dist.family[3] || fam==dist.family[4]) stop(warning("Undifined distribution family"))
  return(CDS)
}

# Treat Normalization
normalize<-function(CDS, bypass, SF) {
  
  if( ! bypass =="bypass"){
    return(estimateSizeFactors(CDS))
  }
  
  pheno$Size_Factor<- SF
  return(CDS)
}


## Select ordering genes 
selectOrderingGenes<-function(CDS, threshold.mean, threshold.fold.dispersion) {
  # Select Genes
  dispTable<-dispersionTable(CDS)
  orderingGenes<-subset(dispTable, mean_expression >= threshold.mean & dispersion_empirical >= threshold.fold.dispersion*dispersion_fit)$gene_id
  # Set CDS
  CDS<-setOrderingFilter(CDS, orderingGenes)
  # Plot genes
  pdf("monocleSelectedGenes.pdf")
  plot_ordering_genes(CDS)
  dev.off()
  # Print selection matrix
  write.table(fData(CDS), col.names=TRUE, row.names=TRUE, sep="\t", file="monocleSelectedGenes.tsv")
  
  return(CDS)
}

### Main script -------------------------------------------------------------------------
## Import data
sce<- import_SCE(args[1])
SF <- BiocGenerics::sizeFactors(sce)

## Make CDS
CDS<-cellDataSet(sce, args[2], as.numeric(args[3]))

## Normalize
CDS<-normalize(CDS, args[4], SF)

# Estimate dispersion
CDS<-estimateDispersions(CDS)

# Select ordering genes
if(as.logical(args[5])) {
  CDS<-selectOrderingGenes(CDS, as.numeric(args[6]), as.numeric(args[7]))
}

# Dimensionality reduction
CDS<-reduceDimension(CDS, max_components=as.numeric(args[8]), reduction_method=args[9], norm_method=args[10])

# Order cells in pseudotime
CDS<-orderCells(CDS, reverse=as.logical(args[11]))

# Plot celltrajectory
pdf("monocleCellsTrajectory.pdf")
color <- pheno[, grep(args[12], colnames(pheno))]
plot_cell_trajectory(CDS, color_by= color)
dev.off()

# Write pheno table
write.table(pData(CDS), file = "monocle_cellMetadata.tsv", col.names=TRUE, row.names=TRUE, sep="\t")

# Save CDS
saveRDS(CDS, "monocle.RDS")
