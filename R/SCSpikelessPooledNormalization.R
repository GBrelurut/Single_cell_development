
#### Load Libraries
library('scran')
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
                       stringsAsFactors=FALSE, sep="\t", row.names=1)
  return(controls)
}


#### Control genes extraction functions =============================================================================================

## Get position of mitochondrial genes in expression matrix row names
get_MT_positions<-function(controls, count){
  ### This function takes as input a control genes table and a count matrix as produced by Eoulsan.
  ### It then seeks features annotated as mitochondrial ("mt_feature") in the count matrix and 
  ### return a vector containing indices of these features.
  
  positions<-c()
  MT_genes<-controls[grep("mt_feature",controls$Type),][,2]
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
  ### It then seeks features annotated as Spike features ("spike_feature") in the count matrix and 
  ### return a vector containing indices of these features.
  
  positions<-c()
  Spike_genes<-controls[grep("spike_feature",controls$Type),][,2]
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

correct_gene_length<-function(controls, count) {
  lengths<-sapply(rownames(count), function(gene){
    controls$Length[grep(gene, rownames(controls))]
  })
  return(count/lengths)
}


#### Plotting functions =====================================================================================================================
  ## Plot PCA
plot.pca<-function(pca, title, color="black"){
  plot(pca$ind$coord[,1:2], 
       pch=19,
       xlab=paste0("PC 1\n(", round(pca$eig[1,2], digits=2), "%)"),
       ylab=paste0("PC 2 (", round(pca$eig[2,2], digits=2), "%)"),
       main=title,
       col=color)
  legend("topright",
         legend=as.character(levels(color)), 
         pch=rep(19, length(levels(color))), 
         col=as.factor(levels(color)), 
         bty="n",
         xpd=NA)
}  
  ## Plot TSNE
plot.tsne<-function(tsne, title, color="black"){
  plot(tsne$Y, 
       pch=19,
       main=title,
       col=color,
       xlab="Dim1",
       ylab="Dim2")
  legend("topright",
         legend=as.character(levels(color)), 
         pch=rep(19, length(levels(color))), 
         col=as.factor(levels(color)), 
         bty="n",
         xpd=NA)
}

  ## Realise all plots for data
plot.data<-function(matrix, color="black", normalized=FALSE) {
  
  # Set Data name
  if(normalized) data<-"Normalized Data"
  else data<-"Raw Data"
  
  # Plot unscaled PCA
  pca<-PCA(t(matrix), graph=FALSE, scale.unit=FALSE)
  pdf(paste0("Scran_", sub(" ", "", data), "_unscaledPCA.pdf"))
  plot.pca(pca, color=color, title=paste0(data, ' PCA (unscaled)'))
  dev.off()
  
  # Plot scaled PCA
  pca<-PCA(t(matrix), graph=FALSE, scale.unit=TRUE)
  pdf(paste0("Scran_", sub(" ", "", data), "_scaledPCA.pdf"))
  plot.pca(pca, color=color, title=paste0(data, ' PCA (scaled)'))
  dev.off()
  rm(pca)
  
  # Plot Tsne on log transformed data
  tsne<-Rtsne(t(log(matrix+1)), pca=FALSE)
  pdf(paste0("Scran_", sub(" ", "", data), "_tsne.pdf"))
  plot.tsne(tsne, color=color, title=paste0(data, ' T-SNE plot \n(log transformed)'))
  dev.off()
}

##### Main function =============================================================================================================================================================================
main<-function(exp_matrix, cell_pheno, genes_table, mode='Nuclear', 
               color.by="Condition", preclustering=FALSE, 
               min_cluster_size=100, force.positive=FALSE,
               cell_cycle=FALSE, organism='mus musculus',
               output1, output2, output3){
  
  ### This function integrates all the parameters given by the user and write the normalized 
  ### expression matrix and the associated size factor in the cells phenotype file.
  
  # Importing Data
  count<-import_exp_matrix(exp_matrix)
  controls<-import_controls(genes_table)
  pheno<-import_pheno(cell_pheno)
  
  # Extracting controls positions 
  MT_positions<-get_MT_positions(controls, count)
  Spike_positions<-get_Spike_positions(controls, count)
  
  # Filter expression matrix
  count<-Filter_exp_matrix(count, MT_positions, Spike_positions, option=mode)
  rm(Spike_positions)
  rm(MT_positions)
    
  # Correct for length if possible and write matrix
  if(! is.null(controls$GeneID)) ccount<-correct_gene_length(controls, count)
  else ccount<-count
  write.table(ccount, paste0(mode, output1), row.names=TRUE, col.names=TRUE, sep='\t')

  # Plotting raw data
  plot.data(ccount, 
            color=as.factor(pheno[,grep(color.by, colnames(pheno))]), 
            normalized=FALSE)
  
  # Constructing SCE object
  SCE<-newSCESet(countData=count)
  
  # Calculate pooling sizes vector
  pooling_sizes<-floor(ncol(count)/20)*10
  pooling_sizes<-seq(from=pooling_sizes/5, to=pooling_sizes, by=pooling_sizes/5)
  
  # Calculating size factors
    if(preclustering){
      clusters<-quickCluster(SCE, min.size=min_cluster_size)
      SCE<-computeSumFactors(SCE, cluster=clusters, sizes=pooling_sizes, positive=force.positive, get.spikes=FALSE, errors=FALSE)
    }else{
      SCE<-computeSumFactors(SCE, cluster=NULL, sizes=pooling_sizes, positive=force.positive, get.spikes=FALSE, errors=FALSE)
    }
  SF<-(sizeFactors(SCE))
  
  # Normalizing data
  SCE<-normalise(SCE, log=FALSE)  
  # treating and saving Cells phenotype file
    # Add size factors
  pheno$Size_factors<-SF
  rm(SF)
    # Annotating for cell cycle
  if(cell_cycle){
    if(organism=='mus musculus'){
      pairs<-readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
    }else{
      if(organism =='homo sapiens'){
        pairs<-readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
      }else{
        stop ("Organism for cell cycle annotation isn't valid.")
      }
    }
    assigned<-cyclone(SCE, pairs) 
    SCE$phase<-rep("S", ncol(SCE))
    SCE$phase[assigned$scores$G1>0.5]<-"G1"
    SCE$phase[assigned$scores$G2M>0.5]<-"G2M"
    SCE$phase[assigned$scores$G1>0.5 & assigned$scores$G2M>0.5]<-"unknown"
    pheno$Cell_cycle<-SCE$phase
  }
  write.table(pheno, output3, row.names=TRUE, col.names=TRUE, sep="\t")

  # Correct for length if possible and save matrix
  if(! is.null(controls$GeneID)) ccount<-correct_gene_length(controls, get_exprs(SCE, "exprs"))
  else ccount<-get_exprs(SCE, "exprs")
  write.table(ccount, paste0(mode, output2), row.names=TRUE, col.names=TRUE, sep='\t')

  
  # Plotting Normalized Data
  plot.data(ccount, 
            color=as.factor(pheno[,grep(color.by, colnames(pheno))]), 
            normalized=TRUE)
}
  
        ### Launch commands
args<-commandArgs(TRUE)
main(args[1], args[2], args[3], args[4], args[5],
      as.logical(args[6]), as.numeric(args[7]), as.logical(args[8]),
      as.logical(args[9]), args[10],
      args[11], args[12], args[13])
