#! /usr/bin/env Rscript

## Import arguments
args<-commandArgs(TRUE)

## Load library
library("methods")
library('scde')

##### Reading and ordering non auxillary tables------------------------------------------

## Import Data
import_SCE <- function(file) {
  return(readRDS(file))
}

## Import data
sce <- import_SCE(args[1])
count <- SummarizedExperiment::assay(sce, "counts")
pheno <- SummarizedExperiment::colData(sce)
group<-factor()

## Treat grouping factor otpions
if(args[3] != "NULL") {
  group<-factor(pheno[,grep(args[3],colnames(pheno))])
} else {
  group<-NULL
}

## Fit error models
count<-apply(count,2,function(x) {storage.mode(x) <- 'integer'; x})
em<-scde.error.models(count, n.cores=as.numeric(args[2]), groups = group, 
                      save.model.plots = as.logical(args[4]),
                      min.nonfailed = as.numeric(args[5]), min.size.entries = as.numeric(args[6]),
                      threshold.segmentation = as.logical(args[7]), min.count.threshold = as.numeric(args[8]),
                      max.pairs = as.numeric(args[9]), min.pairs.per.cell = as.numeric(args[10]),
                      zero.lambda = as.numeric(args[11]), linear.fit = as.logical(args[12]),
                      theta.fit.range=as.numeric(c(args[13], args[14])))

## Write models
write.table(em, args[15], sep="\t", col.names=TRUE, row.names=TRUE)
