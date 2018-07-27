#! /usr/bin/env Rscript

## Import parameters
args<-commandArgs(TRUE)


## load library
library('scde')

## Import Data
import_SCE <- function(file) {
  return(readRDS(file))
}

## Import Model Matrix
import_model_matrix<-function(file){
  ### This function import a matrix, and order it by its row names before returning it.
  ### Name of the file or complete path should be given as input.
  
  em<-read.table(file, header=TRUE, 
                    stringsAsFactors=FALSE, sep='\t', row.names=1) 
  em<-em[order(rownames(em)),]
  return(em)
}

## Import Data
sce <- import_SCE(args[2])
count<- SummarizedExperiment::assay(sce, "counts")
errors<-import_model_matrix(args[1])


## Treat last parameters
quantile<-args[6]
if (quantile == "NULL") { quantile<- NULL
}else {quantile<- as.numeric(quantile)}

max<-args[7]
if (max == "NULL") { max<- NULL
}else {max<-as.numeric(max)}


## Calculate priors according to parameters
if(as.logical(args[4])) png(filename="PriorPlot.png")
prior <- scde.expression.prior(models = errors, counts = count, 
                               length.out = as.numeric(args[3]),
                               show.plot = as.logical(args[4]),
                               pseudo.count = as.numeric(args[5]),
                               max.quantile = quantile,
                               max.value = max )
if(as.logical(args[4])) dev.off()
  
## Write priors
saveRDS(prior, args[8])
