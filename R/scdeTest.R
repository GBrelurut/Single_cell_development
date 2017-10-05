#! /usr/bin/env Rscript


#### sub-Functions =============================================================================

## Import Expression Matrix
import_exp_matrix<-function(file){
  ### This function import a matrix, and order it by its column names before returning it.
  ### Name of the file or complete path should be given as input.
  
  count<-read.table(file, header=TRUE, 
                    stringsAsFactors=FALSE, sep='\t', row.names=1) 
  count<-count[, order(colnames(count))]
  return(count)
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

## Import Prior Matrix
import_prior<-function(file){
  ### This function import a matrix
  prior<-readRDS(file)
  return(prior)
}

## Return Column indices for formula members, members can be experiment specific,
## generic, or a mix of experiment specific and generic
getColumns<-function(formula, pheno, expID){
  expName<-paste0("Exp.", expID,".")
  # Isolate members in a vector
  names<-unlist(strsplit(substr(formula, 2, nchar(formula)), "+", fixed =TRUE))
  # Exclude interaction terms
  names<-names[setdiff(seq(names), grep(":", names))]
  
  # Seek experiment specific columns 
  expCol<-grep(expName, colnames(pheno)) 
  
  # If some columns were found
  if(length(expCol)>0) {
    # Seek experiment specific columns matching required colnames
    namesI<-match(names, sapply(colnames(pheno)[expCol], function(x) {
      return(substr(x, nchar(expName)+1, nchar(x)))
    }), nomatch=0)
    # If all columns are found return indices
    if(length(which(namesI > 0))==length(names)) { return(expCol[namesI])
    # Else seek missing names in other columns
    }else { 
      # If all columns are found return them
      indices<-c(expCol[namesI], (match(names[which(namesI==0)], colnames(pheno), nomatch=0)))
      if(length(which(indices>0)==length(names) ) ) { return (indices)
        # Else stop and print error
      }else { stop(paste0(length(which(indices==0)), " column(s) not found")) }
    }
  }
  # Else seek column in non psecific colnames
  indices<-(match(names, colnames(pheno), nomatch=0))
  # If all columns are found return them
  if(length(which(indices>0)==length(names) ) ) { return (indices)
  # Else stop and print error
  }else { stop(paste0(length(which(indices==0)), " column(s) not found")) }
}

## Return comparison members, in a data frame
getMembers<-function(compa){
  members<-unlist(strsplit(compa, "_vs_", fixed=TRUE))
  if(length(members)>2) stop("Invalid comparisons format")
  submembers1<-unlist(strsplit(members[1], "%", fixed=TRUE))
  submembers2<-unlist(strsplit(members[2], "%", fixed=TRUE))
  if(length(submembers1)!=length(submembers2)) stop("One set is more constrain than the other")
  return(data.frame(set1=submembers1, set2=submembers2))
}

## Return indices for each sets of data, in a list
getLinesIndices<-function(formula, compa, ID, pheno){
  
  # Get indices of columns to use for comparison
  col<-getColumns(formula, pheno, ID)
  
  # Prepare matrix for matching search
  ref<-c()
  if(length(col)>1) ref<-apply(pheno[,col], 2, as.character)
  else ref<-as.character(pheno[,col])
  ref<-rbind(colnames(pheno)[col], ref)
  
  warning(col)
  warning(pheno[,col])
  warning(ref)
          
  # Caculate matching terms
  ref<-apply(ref, 2, function(x) {
    paste0(x[1], x[-1])
  })
  
  #Extract comparison members
  mtable<-getMembers(compa)
  
  warning(mtable$set1)
  warning(mtable$set2)
  
  #Get indices
  
  if(length(col)>1) {
    # For each set
    table<-data.frame(apply(mtable, 2, function(x) {
      # For each line of reference matrix
      logic<-apply(ref, 1, function(y) {
        # For each value in the set : check if the value is present
        sapply(x, function(z){grepl(z, y)})
      })
      # If all value are found in the line, return TRUE, else FALSE
      return(colSums(logic)==length(x))
    }))
  
    #Return indices in a list
    return(list(set1=which(table$set1), set2=which(table$set2)))
  } else {
    table<-data.frame(sapply(mtable, function(x){grepl(x, ref)}))
    return(list(set1=which(table$set1), set2=which(table$set2)))
  }
}

#### Main Script========================================================================

library('scde')

## Import parameters
args<-commandArgs(TRUE)

## Import Data
model<-import_model_matrix(args[1])
count<-import_exp_matrix(args[2])
pheno<-import_pheno(args[3])
prior<-import_prior(args[4])
formula<-args[5]
ID<-args[6]
compa<-args[7]

## Generate comparison factor vector
iList<-getLinesIndices(formula, compa, ID, pheno)
warning(iList$set1)
warning(iList$set2)

groups<-rep(NA, length(nrow(model)))
groups[iList$set1]<-factor("set1")
groups[iList$set2]<-factor("set2")
warning(groups)

## Get batch vector
if(args[10] != "NULL") { batchCol<-as.factor(pheno[, grep(args[10], colnames(pheno))])
}else { batchCol<-NULL }

## Launch test
result<-scde.expression.difference(model, count, prior, groups, n.cores=as.numeric(args[8]),
                           n.randomizations=as.numeric(args[9]), batch=batchCol,
                           return.posteriors = as.logical(args[11]))



## Print results
if(as.logical(args[11])) {
  diff.post<-result$difference.posterior
  joint.post<-result$joint.posteriors
  results<-result$results
  
  write.table(joint.post[[1]], paste0("Exp-", ID, "-JointPosteriors",
                                      names(joint.post)[1],".tsv"),
              col.names=TRUE, row.names=TRUE, sep="\t")
  write.table(joint.post[[2]], paste0("Exp-", ID, "-JointPosteriors",
                                      names(joint.post)[2],".tsv"),
              col.names=TRUE, row.names=TRUE, sep="\t")
  write.table(diff.post, paste0("Exp-", ID, "-DifferencialExpressionPosteriors.tsv"), 
              col.names=TRUE, row.names=TRUE, sep="\t")
  
  rm(joint.post)
  rm(diff.post)
  gc()
  
} else {
  results<-result$results
}

# Return batch correction
if(!is.null(batchCol)) {
  batch.biased<-result$results
  batch.effect<-result$batch.effect
  write.table(batch.effect, paste0("Exp-", ID, "-BatchEffect.tsv"), col.names=TRUE, row.names=TRUE, sep="\t")
  write.table(batch.biased, paste0("Exp-", ID, "-BatchUncorrectedResult.tsv"), col.names=TRUE, row.names=TRUE, sep="\t")
  
  rm(batch.biased)
  rm(batch.effect)
  
  if(as.logical(args[11])) {
    batch.diff<-result$batch.adjusted.difference.posteriors
    write.table(batch.diff, paste0("Exp-", ID, "BatchAdjustedPosteriors.tsv"))
  }
  gc()
  
  results<-result$batch.adjusted
  
}


results$cProbability<-pnorm(abs(results$cZ), lower.tail=FALSE)*2

write.table(results, file=args[12], sep="\t", row.names=TRUE, col.names=TRUE)
