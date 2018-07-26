#! /usr/bin/env Rscript

library('ggplot2')
library('SingleCellExperiment')

##### Data pre-processing ===============================================================

##### Reading and ordering non auxillary tables------------------------------------------

## Import Expression Matrix
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
  spikes <- rowData(SCE)$Type == "spike"
  SCE <- SCE[! spike,]
  
  if(mode == "Endogenous") return(SCE)
  
  # Get and remove mitochondrial genes
  mt <- rowData(SCE)$Type == "mitochondrial"
  return(SCE[! mt,])
  
}

#### Calcul and filter functions======================================================================

#### Get metrics-----------------------------------------------------------------------------------

## Calculate proportion of mapped reads
### This function takes as input a single cell experiment and a type of features.
### It then calculate the proportion of reads mapped to the corresponding positions.
get_prop <- function(sce, type = "mitochondrial") {
  total <- colSums(SummarizedExperiment::assay(sce, "counts"))
  ctrls <- colSums(SummarizedExperiment::assay(sce[colData(sce)$Type == type,], "counts"))
  return(ctrls/total)
}

#### Filters------------------------------------------------------------------------------

## Filter cell by number of failure
filter_cells<-function(cells, threshold){
  ### This function takes as input a vector with repeated element and an integer value for threshold.
  ### It then detects and returns elements appearing more than the given threshold.
  
  cells<-data.frame(table(cells))
  filt<-which(cells$Freq>=threshold)
  return(as.numeric(levels(cells[,1])[filt]))
}

#### Main function===============================================================================

main<-function(file1,
               detection=10, exp_option='Nuclear',
               exp_threshold=4000, nbr_threshold=500000,
               propmt_threshold=0.2,
               propsp_threshold=0.5,
               nb_filters=1,
               output1, output2){
  ### This function integrates all parameters given by user and realises the quality 
  ### checking of the data.
  ### Set threshold to 0 to disable filtering
  
  # Importing Data-------------------------------------------  
  sce <- import_SCE(file1)

  
  # Initializing filtering result vectors--------------------
  mt_cells<-c()
  sp_cells<-c()
  
  # Treating data for advanced metrics ----------------------
  if(propmt_threshold){
    prop <- get_prop(sce, "mitochondrial")
    SummarizedExperiment::colData(sce)$Prop_Mt <- prop
    mt_cells <- which(prop > propmt_threshold)
  }
  
  if(propsp_threshold){
    prop <- get_prop(sce, "spike")
    SummarizedExperiment::colData(sce)$Prop_Sp <- prop
    sp_cells <- which(prop > propsp_threshold)
  }
  
  # Treating data for basic metrics -------------------------
  fmatrix <- SingleCellExperiment::assay(filter_SCE(sce, mode = exp_option), "counts")
  exp_features<-colSums(fmatrix > detection)
  nb_reads<-colSums(fmatrix)
  rm(fmatrix)
  
  if(exp_threshold){
    exp_cells <- which(exp_features < exp_threshold)
  }
  if(nbr_threshold){
    nbr_cells<-which(nb_reads< nbr_threshold)
  }
  SummarizedExperiment::colData(sce)$Nb_features<-exp_features
  rm(exp_features)
  SummarizedExperiment::colData(sce)$Nb_reads<-nb_reads
  rm(nb_reads)
  
  # Code from previous version allowing different options for feature and reads counting------  
  # fmatrix<-Filter_exp_matrix(count, MT_positions, Spike_positions, option=exp_option)
  # exp_features<-get_detected_genes(fmatrix, detection)
  # rm(fmatrix)
  # if(exp_threshold){
  #   exp_cells<-which(exp_features < exp_threshold)
  # }
  # pheno$Nb_features<-exp_features
  # rm(exp_features)
  # fmatrix<-Filter_exp_matrix(count, MT_positions, Spike_positions, option=nbr_option)
  # nb_reads<-colSums(fmatrix)
  # rm(fmatrix)
  # if(nb_reads){
  #   nbr_cells<-which(nb_reads< nbr_threshold)
  # }
  # pheno$Nb_reads<-nb_reads
  # rm(nb_reads)
  
  
  # Export results ------------------------------------------
  saveRDS(sce, output1)
  
  # Extract filtered cells -----------------------------------
  filter<-filter_cells(c(exp_cells, nbr_cells, mt_cells, sp_cells), nb_filters)
  
  
  
  # Plot Raw data with filtered cells annotated -------------- 
  pdf(paste0('Raw_Cellplot.pdf')) ##### add function to get experiment from input files
  plot(SummarizedExperiment::colData(sce)$Nb_reads, 
       SummarizedExperiment::colData(sce)$Nb_features, xlab='Number of reads', ylab='Number of detected features', 
       xlim=c(min(pheno$Nb_reads)-500,max(pheno$Nb_reads)+1000),
       ylim=c(min(pheno$Nb_features)-50, max(pheno$Nb_features)),
       pch=19, col=alpha("gray", 0.5))
  points(x=pheno$Nb_reads[filter], y=pheno$Nb_features[filter], col='red', pch=19)
  dev.off()
  
  # Filtered and export table and matrices--------------------
  if (length(filter) > 0L) {
    sce <- sce[,-filter]
  }
  saveRDS(sce, output2)
  
  # Export filtered dot plot ---------------------------------
  pdf(paste0( 'Filtered_cellplot.pdf')) ##### add function to get experiment from input files
  plot(SummarizedExperiment::colData(sce)$Nb_reads, 
       SummarizedExperiment::colData(sce)$Nb_features, 
       xlab='Number of reads', ylab='Number of detected features',
       pch=19, col=alpha("gray", 0.5))
  dev.off()
}

args<-commandArgs(TRUE)

main(args[1],
     detection=as.numeric(args[2]),
     exp_option = args[3],
     exp_threshold = as.numeric(args[4]),
     nbr_threshold = as.numeric(args[5]),
     propmt_threshold = as.numeric(args[6]), 
     propsp_threshold = as.numeric(args[7]), 
     nb_filters = as.numeric(args[8]),
     args[9], args[10])
