#! /usr/bin/env Rscript

## Import arguments
args<-commandArgs(TRUE)

## Load library
library('monocle')

### Functions ---------------------------------------------------------------------------
## Import Expression Matrix
import_exp_matrix<-function(file){
  ### This function import a matrix, and order it by its column names before returning it.
  ### Name of the file or complete path should be given as input.
  
  count<-read.table(file, header=TRUE, 
                    stringsAsFactors=FALSE, sep='\t', row.names=1) 
  count<-count[, order(colnames(count))]
  return(count)
}


##Import Reduced Design File
import_pheno<-function(file){
  ### This function import a table and named its row by its second column before ordering them and returning it.
  ### Name of the file or complete path should be given as input.
  
  pheno<-read.table(file, header=TRUE, 
                    stringsAsFactors=FALSE, sep='\t') 
  rownames(pheno)<-pheno[,2]
  pheno<-pheno[ order(rownames(pheno)),]
  return(pheno)
}

## Create Cell Data Set
cellDataSet<-function(count, fam, detection.threshold) {
  ## Set option verctor
  dist.family<-c("negbinomial.size", "negbinomial", "tobit", "gaussianff") # Last families are untested, but conserved to help module update
  # Check parameter
  if(! fam %in% dist.family) stop(warning("Undifined distribution family"))
  # Create CDS
  if(fam==dist.family[1]) CDS<-newCellDataSet(as.matrix(count), expressionFamily=negbinomial.size(), lowerDetectionLimit=detection.threshold)
  if(fam==dist.family[2]) CDS<-newCellDataSet(as.matrix(count), expressionFamily=negbinomial(), lowerDetectionLimit=detection.threshold)
  if(fam==dist.family[3] || fam==dist.family[4]) stop(warning("Undifined distribution family"))
  return(CDS)
}

# Treat Normalization
normalize<-function(CDS, bypass, pheno) {
  
  if( ! bypass=="bypass"){
    pData(CDS)<-pheno
    return(estimateSizeFactors(CDS))
  }
  
  pheno$Size_Factor<-pheno$Size_factors
  pData(CDS)<-pheno
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
count<-import_exp_matrix(args[1])
pheno<-import_pheno(args[2])

## Make CDS
CDS<-cellDataSet(count, args[3], as.numeric(args[4]))

## Normalize
CDS<-normalize(CDS, args[5], pheno)

# Estimate dispersion
CDS<-estimateDispersions(CDS)

# Select ordering genes
if(as.logical(args[6])) {
  CDS<-selectOrderingGenes(CDS, as.numeric(args[7]), as.numeric(args[8]))
}

# Dimensionality reduction
CDS<-reduceDimension(CDS, max_components=as.numeric(args[9]), reduction_method=args[10], norm_method=args[11])

# Order cells in pseudotime
CDS<-orderCells(CDS, reverse=as.logical(args[12]))

# Plot celltrajectory
pdf("monocleCellsTrajectory.pdf")
plot_cell_trajectory(CDS, color_by= as.character(args[13]))
dev.off()

# Write pheno table
write.table(pData(CDS), file=args[2], col.names=TRUE, row.names=TRUE, sep="\t")
