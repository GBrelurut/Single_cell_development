#! /usr/bin/env Rscript

library('scales')
library('parallel')
library('SingleCellExperiment')

#### Data pre-processing ===============================================================

#### Reading and ordering non auxillary tables------------------------------------------

## Import Data
## Import SCE
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


#### Saturation functions ========================================================================
## Resample function
resample<-function(vector, iterations=10, threshold=10, round=10000, nbPoints=20) {
  ### This function takes as input a vector of count for various element (one element in 
  ### the vector = count for one element)
  ### Options are: 
  ###  iterations:  number of iteration for each step, default value is 10
  ###  threshold:   number of event for an element to be considered as detected, 
  ###               default is 10
  ###  round :      Magnitude for down rounding counts, default is 10 000
  ###  nbPoints :   Number of points to sample
  ### returns : sampling results in a dataframe
  
  # getting down rounded max depth of resampling
  tot_count<-sum(vector)
  tot_count<-floor(tot_count/round)*round
  
  # Creating vector for random sampling
  data<-rep(1:length(vector), vector)
  
  # Set increment
  steps=tot_count/nbPoints
    
  # Create sampling vector
  depths<-seq(0, tot_count, steps)
  
  # Internal function
  ## Draw and count (sampling)
  drawNcount<-function(index, data, depth, threshold){
    samp<-data.frame(table(sample(data, depth, replace=FALSE)))
    count<-sum(samp$Freq>=threshold)
    names(count)<-paste0('draw', index)
    return(count)
  }
  ## Resume sampling
  resume<-function(depth, data, threshold, iterations){
    sampl<-sapply(1:iterations, drawNcount, data=data, 
                  depth=depth, threshold=threshold)
    return(c(depth, mean(sampl), sd(sampl)))
  }
  
  # Resample
  result<-sapply(depths, resume, data=data, threshold=threshold, iterations=10)
  
  #return result
  rownames(result)<-c('depths','means','sd')
  return(data.frame(t(result)))
}


## Fitting function
fitMMModel<-function(table){
  ### This function fits a Michaelis Menten model using a table of resampled features
  ### table :   table of resampled features
  ### return :  a List : max = maximum of the Model,  half = X value to have 1/2 of max
  ###           min = minimum of the fitting function
  
  # function to minimize for non linear model
  michaelis<-function(param, data){
    max<-param[1]
    half<-param[2]
    diff<-sum(abs(data$means-(max*data$depths/(data$depths+half))))
    return(diff)
  }
  # Fit a linear model using Hanes-Woolf equation to obtain near optimal parameters value
  HW<-lm(means~depths, data.frame(depths=table[-which(table[,2]==0),1], 
                                  means=table[-which(table[,2]==0),1]/table[-which(table[,2]==0),2]))
  
  # Extract parameters value
  HW.max<-1/HW$coefficients[2]
  HW.half<-HW$coefficients[1]/HW$coefficients[2]
  
  # Fit non linear model using near optimal starting values
  results<-nlm(michaelis, c(HW.max, HW.half), table)
  return(list(max=results$estimate[1], half=results$estimate[2], min=results$minimum) )
}


## Fucntion for fitting evaluation
testFit<-function(measured, predicted, name){
  ### Fit and plot a linear model such as predicted = f(measured)
  ### measured :  Measured values (numeric vector)
  ### predicted : values predicted by model (numeric vector)
  ### return : R² of linear fit (indicating goodness of fit)
  
  # Fit linear model
  model<-lm(predicted~measured)
  # Get R²
  r<-summary(model)$r.squared
  
  # Plot and save 
  pdf(paste0(name, 'FittingPlot.pdf'))
  plot(measured, predicted, xlab='Measured values', ylab='Predicted values', 
       main=paste0(name, ' goodness of fit'), pch=19)
  abline(model, col='chocolate3')
  dev.off()
  
  #Return R²
  return(r)
}


## Plotting function
plotSaturation<-function(table, predicted, name) {
  ### This function plots and saves saturation resampling given in a dataframe
  ### Parameters : 
  ###     table : table gathering data to plot with sampling on first column means 
  ###             on second column and sd on third column
  ###     name : name of the sample (for plot title)
  
  sampl<-table[,1]
  mean<-table[,2]
  sd<-table[,3]
  
  # data
  plot<-ggplot(data.frame(reads=sampl, genes=mean), 
               aes(x=reads, y=genes))
    # Points and line
    plot<-plot+geom_line(colour='royalblue3')
    plot<-plot+geom_point(colour='royalblue1',size=2)
    # Error bars
    plot<-plot+geom_errorbar(aes(ymin=mean-sd, 
                      ymax=mean+sd), colour='royalblue3', width=(sampl[2]-sampl[1])/5)
    # Adding model curve
    plot<-plot+geom_line(aes(x=sampl, y=predicted), colour='chocolate3')
    # esthetics
    plot<-plot+theme_bw()
    plot<-plot+xlab('Number of reads')+ylab('Number of detected features')
    plot<-plot+ggtitle(paste0(name,' Saturation Curve '))
  ggsave(paste0('Saturation_curve_', name, '.pdf'), plot=plot, device='pdf')
}


## Function to treat columns from modified counting matrix
saturationCol<-function(col, iterations=10, threshold=10, round=10000, nPoints=20) {
  ###This function takes as input an expressionvector and gives as output a
  ###saturation plot, fitting evaluation plot, and a vector containing results 
  ###of the fitting. Requires ggplot2.
  ###Options are: 
  ### iterations: number of iteration for each step, default value is 10
  ### threshold:  number of event for a gene to be considered as detected
  ###             (set as a parameter for the module, default is 10)
  ### round:      Magnitude for down rounding counts data, default is 
  ###             10 000.
  ### nPoints :   number of points to sample for saturation
  ### Return :    a vector containing : model coefficients (maximum and X  
  ###             value for Y=1/2*maximum), fitting value, saturation estimate
  
  # Extract count and name
  l<-length(col)
  name<-col[l]
  exp<-as.numeric(col[-l])
  
  # Resample data
  rs<-resample(exp, threshold=threshold, iterations=iterations, round=round, nbPoints=nPoints)
  # Check sampling, if cell has not enough reads return NA values
  if(nrow(rs)< nPoints) return(rep(NA, 4))
  
  #fit model on resampling
  MM<-fitMMModel(rs)
  # Get predicted values
  predicted<-MM$max*rs$depths/(MM$half+rs$depths)
  
  # Plotdata
  if(! dir.exists('SaturationPlots')) dir.create('SaturationPlots')
  setwd('SaturationPlots')
  plotSaturation(rs, predicted, name)
  setwd('..')
  
  # Test Goodness of fit
  if(! dir.exists('FitEvaluation')) dir.create('FitEvaluation')
  setwd('FitEvaluation')
  r<-testFit(rs$means, predicted, name)
  setwd('..')
  
  # Evaluate saturation
  sat<-sum(exp >= threshold)/MM$max
  
  # Return results
  return(c(MM$max, MM$half, sat, r))
}

## Function for index correction
correctIndex<-function(index, ref){
  ### given an index and removed indices, correct the index
  ### index : a numeric value
  ### ref :   a numeric vector containing removed values (increasing)
  ### return corrected index value
  i<-1
  while(ref[i] <= index & i <=length(ref)){
    i<-i+1
    index<-index+1
  }
  return(index)
}

## Function identifying outliers (won't work if more than 5000 values are passed)
getOutliers<-function(values, normThreshold=0.001){
  ### This function identify long left tail outliers
  ### values : numerical vector
  ### normThreshold : probability under gaussian expectation to reach
  ### return : vector containing indices of removed values
  
  # Initialize variables
  p<-shapiro.test(values)$p.value
  proba<-p
  rmv<-NULL
  
  # Remove outliers
  while(p<normThreshold) {
    if( ! is.null(rmv)) {
      mini<-which.min(values[-rmv])
      mini<-correctIndex(mini, rmv)
      rmv<-c(rmv, mini)
      rmv<-rmv[order(rmv, decreasing=FALSE)]
    }
    else {
      mini<-which.min(values)
      rmv<-mini
    }
    if( max(values[-rmv]) != min(values[-rmv]) & length(rmv)<length(values)-3) {
      p<-shapiro.test(values[-rmv])$p.value
      proba<-c(proba, p)
    } else {
      warning("Outliers removal failed")
      break
    }
  }
  
  
  # Plot P values
  pdf('gaussianProbability.pdf')
  plot(1:length(proba)-1, proba, type='b', pch=19, 
       xlab='Number of removal', ylab='P value (Shapiro)')
  dev.off()
  
  # return indices
  return(rmv)
}


## Saturation global function
saturation <- function (count, satThreshold=0.7, fitThreshold=0.97, iterations = 10, 
                        threshold=10, round=10000, nPoints=20, nCores=2) {
  ###This function takes as input an expression matrixand gives as output a
  ###saturation plot and a fitting plot for heach column of the matrix (Requires ggplot2)
  ###
  ###Options are: 
  ### iterations: number of iteration for each step, default value is 10
  ### threshold:  number of event for a gene to be considered as detected
  ###             (set as a parameter for the module, default is 10)
  ### round:      Magnitude for down rounding counts data, default is 
  ###             10 000.
  ### nPoints :   number of points to sample for saturation
  ### nCores :    number of threads to run
  ### Return :    a vector containing : model coefficients (maximum and X  
  ###             value for Y=1/2*maximum), fitting value, saturation estimate
  
  # Add a line containing Sample name to input matrix
  mat<-rbind(count, colnames(count))
  
  # Prepare parallelized computing
  clus<-makeCluster(nCores)
  clusterEvalQ(clus, library('ggplot2'))
  clusterExport(clus, c('iterations', 'threshold', 'round', 'nPoints'), envir=environment())
  clusterExport(clus, c( 'resample', 'fitMMModel', 'plotSaturation', 'testFit', 'saturationCol'))
  
  # Following operations are applied to each cell
  metrics<-parCapply(clus, mat, saturationCol, iterations=iterations, threshold=threshold,
                     round=round, nPoints=nPoints )
  stopCluster(clus)
  metrics<-matrix(metrics, nrow=4, byrow=FALSE)
  colnames(metrics)<-colnames(count)
  rownames(metrics)<-c('max', 'half', 'saturation', 'r.squared')
  
  # Change Working Directory to save plots
  setwd('FitEvaluation')
  
  # Get cells to suppress indices
  i<-which(apply(metrics[1:2,], 2, function(x){any(x <0, na.rm=TRUE)}) | metrics[3,] < satThreshold | apply(metrics, 2, anyNA))
  outliers<-which(metrics[4,] < fitThreshold) 
  i<-c(i, outliers)
  
  # Plot data before and after removal
  rmv<-rep('grey', ncol(metrics))
  rmv[i]<-'red'
  
    ## Full set
  pdf('qualityAll.pdf')
  plot(t(metrics[c(1,4),]), xlab='Maximum', ylab='R squared', pch=19, 
       col=rmv, main='Saturation quality \n (all cells)')
  dev.off()
  
  pdf('densityAll.pdf')
  plot(density(metrics[4,], na.rm=TRUE), main="R squared density \n (all cells)")
  dev.off()
  
    ## Filtered cell
  if(length(i) > 0) {
    pdf('qualityFiltered.pdf')
    plot(t(metrics[c(1,4),-i]), xlab='Maximum', ylab='R squared', pch=19, main='Saturation quality \n (good cells)')
    dev.off()
  
    pdf('densityFiltered.pdf')
    plot(density(metrics[4, -i], na.rm=TRUE), main="R squared density \n (good cells)")
    dev.off()
  }
  
  # Back to main directory
  setwd('..')
  
  # return indices  
  return(list(indices=i, metrics=t(metrics)))
}


#### Main function===============================================================================
main<-function(file, 
               detection=10,  exp_option='Nuclear',
               satThreshold=0.7, fitThreshold=0.01, round=10000, nCores=2,
               propmt_threshold=0.2, propsp_threshold=0.5, nb_filters=1,
               output1, output2){
  ### This function integrates all parameters given by user and realises the quality 
  ### checking of the data.
  ### Set threshold to 0 to disable filtering
  
  # Importing Data-------------------------------------------  
  sce <- import_SCE(file)
  
  # Initializing filtering result vectors--------------------
  mt_cells<-c()
  sp_cells<-c()
  
  # Treating data for advanced metrics ----------------------
  if(propmt_threshold){
    prop <- get_prop(sce, "mitochondrial")
    rowData(sce)$Prop_Mt <- prop
    mt_cells <- which(prop > propmt_threshold)
  }
  
  if(propsp_threshold){
    prop <- get_prop(sce, "spike")
    rowData(sce)$Prop_Sp <- prop
    sp_cells <- which(prop > propsp_threshold)
  }
  
  # Filtering SCE ----------------------------------------
  fmatrix <- SingleCellExperimen::assay(filter_SCE(sce, mode = exp_option), "counts")
  
  # Extracting basic metrics --------------------------------
  exp_features<- colSums(fmatrix > detection)
  nb_reads <- colSums(fmatrix)
  
  # Saving in phenotype table -------------------------------
  SummarizedExperiement::colData(sce)$Nb_features <- exp_features
  rm(exp_features)
  SummarizedExperiment::colData(sce)$Nb_reads<-nb_reads
  rm(nb_reads)
   
  # Finding unsaturated cells -------------------------------
  unsat<-saturation(fmatrix, threshold = detection, satThreshold = satThreshold, 
                      fitThreshold = fitThreshold, round = round, nCores= nCores)
  
  colData(sce) <- cbind(colData(sce), unsat$metrics)
  
  # Export results ------------------------------------------
  saveRDS(sce, output1)
  
  # Extract filtered cells -----------------------------------
  filter<-filter_cells(c(unsat$indices, mt_cells, sp_cells), nb_filters)
 
  # Plot Raw data with filtered cells annotated --------------
  colour<-rep(alpha('gray', 0.5), nrow(SummarizedExperiment::colData(sce)))
  colour[filter]<-'red'
  pdf(paste0('Raw_Cellplot.pdf'))
  plot(SummarizedExperiment::colData(sce)$Nb_reads, 
       SummarizedExperiment::colData(sce)$Nb_features, 
       xlab='Number of reads', ylab='Number of detected features', 
       xlim=c(min(pheno$Nb_reads)-1000,max(pheno$Nb_reads)+1000),
       ylim=c(min(pheno$Nb_features)-50, max(pheno$Nb_features)+50),
       pch=19, col=colour)
  dev.off()
  
  # Filtered and export table and matrices--------------------
  if (length(filter) > 0) {
    sce < sce[ , - filter]
  }
  saveRDS(sce, output2)
  
  # Export filtered dot plot ---------------------------------
  pdf(paste0( 'Filtered_cellplot.pdf')) 
  plot(SummarizedExperiment::colData(sce)$Nb_reads, 
       SummarizedExperiment::colData(sce)$Nb_features, 
       xlab='Number of reads', ylab='Number of detected features',
       pch=19, col=alpha('gray', 0.5))
  dev.off()
}

#### Running script==============================================================================
args<-commandArgs(TRUE)

main(args[1],  
     detection=as.numeric(args[2]), exp_option = args[3],
     satThreshold = as.numeric(args[4]), fitThreshold = as.numeric(args[5]),
     round = as.numeric(args[6]), nCores = as.numeric(args[7]),
     propmt_threshold = as.numeric(args[8]), 
     propsp_threshold = as.numeric(args[9]), 
     nb_filters = as.numeric(args[10]),
     args[11], args[12])
