
#### Load Libraries
library('FactoMineR')
library('Rtsne')

##### Data pre-processing ==================================================================================================
  ##### Reading and ordering non auxillary tables

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
                    stringsAsFactors=FALSE, sep='\t', row.names=1) 
  pheno<-pheno[ order(rownames(pheno)),]
  return(pheno)
}

        ## Import control genes table
import_controls<-function(file){
  controls<-read.table(file, header=TRUE,
                       stringsAsFactors=FALSE, sep='\t',
                       row.names=1)
  return(controls)
}


#### Control genes extraction functions =============================================================================================

## Get position of mitochondrial genes in expression matrix row names
get_MT_positions<-function(controls, count){
  ### This function takes as input a control genes table and a count matrix as produced by Eoulsan.
  ### It then seeks features annotated as mitochondrial ('mt_feature') in the count matrix and 
  ### return a vector containing indices of these features.
  
  positions<-c()
  MT_genes<-controls[grep('mt_feature',controls$Type),][,2]
  if (length(MT_genes>0)){
    for (gene in 1:length(MT_genes)){
      name<-MT_genes[gene]
      positions<-c(positions, grep(name, rownames(count)))
    }
  }
  return(as.numeric(positions))
}

## Get position of spike-ins genes in expression matrix row names
get_Spike_positions<-function(controls, count){
  ### This function takes as input a control genes table and a count matrix as produced by Eoulsan.
  ### It then seeks features annotated as Spike features ('spike_feature') in the count matrix and 
  ### return a vector containing indices of these features.
  
  positions<-c()
  Spike_genes<-controls[grep('spike_feature',controls$Type),][,2]
  if (length(Spike_genes>0)){
    for (gene in 1:length(Spike_genes)){
      name<-Spike_genes[gene]
      positions<-c(positions, grep(name, rownames(count)))
    }
  }
  return(as.numeric(positions))
}

  #### Processing Expression Matrix

        ## Expression Matrix filtering
Filter_exp_matrix<-function(count, Mt_positions, Spike_positions, option='Nuclear'){
  ### This function takes as input a count matrix and positions to be treated.
  ### It then process the matrix according to three options (see below), removing 
  ### unwanted genes and returning the filtered matrix.
  ### Options : 
  ###           -Endogenous : Matrix is filtered for exogenous genes (spike-ins)
  ###           - Nuclear :   Matrix is filtered for exogenous genes and mitochondrial genes
  
  
    if(grepl('Nuclear', option, fixed=TRUE)){
      if(length(Mt_positions)>0 || length(Spike_positions)>0){
        count<-count[-(c(Mt_positions, Spike_positions)),]
      }else{
        warning('Positions if of length O. \n Nothing was filtered')
      }
    }else{
      if(grepl('Endogenous', option, fixed=TRUE)){
        if(length(Spike_positions)>0){
          count<-count[-c(Spike_positions),]
        }
      }else{
        stop('Invalid option')
      }
    }
  return(count)
}

## Function for length correction
correct_gene_length<-function(controls, count) {
  lengths<-sapply(rownames(count), function(gene){
    controls$Length[grep(gene, rownames(controls))]
  })
  return(count/lengths)
}


## SF function
calculate_SF<-function(count, mode='CPM', genes=NULL ) {
  if(mode == 'TPM' && ! is.null(genes)) count <- correct_gene_length(genes, count)
  return(colSums(count)/1e6)
}

#### Plotting functions =====================================================================================================================
  ## Plot PCA
plot.pca<-function(pca, title, color='black'){
  plot(pca$ind$coord[,1:2], 
       pch=19,
       xlab=paste0('PC 1\n(', round(pca$eig[1,2], digits=2), '%)'),
       ylab=paste0('PC 2 (', round(pca$eig[2,2], digits=2), '%)'),
       main=title,
       col=color)
  legend('topright',
         legend=as.character(levels(color)), 
         pch=rep(19, length(levels(color))), 
         col=as.factor(levels(color)), 
         bty='n',
         xpd=NA)
}  
  ## Plot TSNE
plot.tsne<-function(tsne, title, color='black'){
  plot(tsne$Y, 
       pch=19,
       main=title,
       col=color,
       xlab='Dim1',
       ylab='Dim2')
  legend('topright',
         legend=as.character(levels(color)), 
         pch=rep(19, length(levels(color))), 
         col=as.factor(levels(color)), 
         bty='n',
         xpd=NA)
}

  ## Realise all plots for data
plot.data<-function(matrix, color='black', normalized=FALSE) {
  
  # Set Data name
  if(normalized) data<-'Normalized Data'
  else data<-'Raw Data'
  
  # Plot unscaled PCA
  pca<-PCA(t(matrix), graph=FALSE, scale.unit=FALSE)
  pdf(paste0('Sum_', sub(' ', '', data), '_unscaledPCA.pdf'))
  plot.pca(pca, color=color, title=paste0(data, ' PCA (unscaled)'))
  dev.off()
  
  # Plot scaled PCA
  pca<-PCA(t(matrix), graph=FALSE, scale.unit=TRUE)
  pdf(paste0('Sum_', sub(' ', '', data), '_scaledPCA.pdf'))
  plot.pca(pca, color=color, title=paste0(data, ' PCA (scaled)'))
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
  pdf(file=paste0('Sum_', sub(' ', '', data), '_tsne.pdf'), width=8*length(n.neighbors), height=8)
  layout(matrix(c(1:length(n.neighbors)), nrow=1, byrow=TRUE))
  for(nn in n.neighbors) {
    tsne<-Rtsne(t(log(matrix+1)), pca=FALSE, perplexity=nn)
    plot.tsne(tsne, color=color, title=paste0(data, ' T-SNE plot \n(log transformed, perplexity = ', nn, ')'))
  }
  dev.off()
  layout(matrix(c(1), nrow=1, byrow=TRUE))
}

##### Main function =============================================================================================================================================================================
main<-function(exp_matrix, cell_pheno, genes_table, mode='Nuclear', 
               color.by='Condition', norm='CPM', 
               length_correction=TRUE,
               output){
  
  ### This function integrates all the parameters given by the user and write the normalized 
  ### expression matrix and the associated size factor in the cells phenotype file.
  
  # Importing Data
  count<-import_exp_matrix(exp_matrix)
  controls<-import_controls(genes_table)
  pheno<-import_pheno(cell_pheno)
  
  # Extracting controls positions 
  MT_positions<-get_MT_positions(controls, count)
  Spike_positions<-get_Spike_positions(controls, count)
  
  # Treating expression matrix
  count<-Filter_exp_matrix(count, MT_positions, Spike_positions, option=mode)
  rm(Spike_positions)
  rm(MT_positions)

  # Correct for length if possible and save expression matrix
  if(length_correction && ! is.null(controls$GeneID)){
    ccount<-correct_gene_length(controls, count)
    write.table(ccount,  'lengthCorrectedCount.tsv', row.names=TRUE, col.names=TRUE, sep='\t')
  }
  else {
    ccount<-count
  } 
  
  
  # Plotting raw data
  plot.data(count, 
            color=as.factor(pheno[,grep(color.by, colnames(pheno))]), 
            normalized=FALSE)
  

  # Calculating size factors
  genes <- NULL
  if(norm == 'TPM') {
    if(! is.null(controls$GeneID)) {
      genes <- controls
    } else {
      warning('TPM normalization is not available for gene counting, shifting to CPM')
      norm <- 'CPM'
    }
  }
  SF<-calculate_SF(count, norm, genes)
  
  # Normalizing data
  count<-t(t(count)/SF)
  
  # treating and saving Cells phenotype file
  pheno$Size_factors<-SF
  rm(SF)
  write.table(pheno, output, row.names=TRUE, col.names=TRUE, sep='\t')

  # Correct for length if possible and write table
  if(length_correction && ! is.null(controls$GeneID)){
    ccount<-correct_gene_length(controls, count)
    write.table(ccount,  'lengthCorrectedNormalizedCount.tsv', row.names=TRUE, col.names=TRUE, sep='\t')
  }
  else {
    ccount<-count
    write.table(ccount,  'normalizedCount.tsv', row.names=TRUE, col.names=TRUE, sep='\t')
  } 
  
  # Plotting Normalized Data
  plot.data(ccount, 
            color=as.factor(pheno[,grep(color.by, colnames(pheno))]), 
            normalized=TRUE)
  
  # Clean directory
  if(file.exists('Rplots.pdf')) file.remove('Rplots.pdf')
  
}
  
        ### Launch commands
args<-commandArgs(TRUE)
main(args[1], args[2], args[3], args[4], args[5],
      args[6], as.logical(args[7]),
      args[8])
