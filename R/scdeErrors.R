#! /usr/bin/env Rscript

## Import arguments
args<-commandArgs(TRUE)

## Load library
library('scde')

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
                    stringsAsFactors=FALSE, sep='') 
  rownames(pheno)<-pheno[,2]
  pheno<-pheno[ order(rownames(pheno)),]
  return(pheno)
}

## Import data
count<-import_exp_matrix(args[1])
pheno<-import_pheno(args[2])
group<-factor()

## Treat grouping factor otpions
if(args[4] != "NULL") {
  group<-factor(pheno[,grep(args[4],colnames(pheno))])
} else {
  group<-NULL
}

## Fit error models
em<-scde.error.models(count, n.cores=as.numeric(args[3]), groups = group, 
                      save.model.plots = as.logical(args[5]),
                      min.nonfailed = as.numeric(args[6]), min.size.entries = as.numeric(args[7]),
                      threshold.segmentation = as.logical(args[8]), min.count.threshold = as.numeric(args[9]),
                      max.pairs = as.numeric(args[10]), min.pairs.per.cell = as.numeric(args[11]),
                      zero.lambda = as.numeric(args[12]), linear.fit = as.logical(args[13]),
                      theta.fit.range=as.numeric(c(args[14], args[15])))

## Write models
write.table(em, args[16], sep="\t", col.names=TRUE, row.names=TRUE)
