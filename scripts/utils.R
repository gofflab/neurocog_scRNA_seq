library(cellrangerRkit)
library(stringr)
library(monocle)
library(reticulate)
use_condaenv('r-reticulate')
import("louvain")
library(VGAM)
library(pheatmap)
library(gridExtra)
library(dplyr)
library(NNLM)
library(viridis)
library(reshape2)
library(DT)

# Below must be run once, and then is not needed again
#conda_install('r-reticulate','umap-learn')#, pip = T, pip_ignore_installed = T)
#conda_install('r-reticulate','louvain')#, pip = T, pip_ignore_installed = T)
#devtools::install_github('VPetukhov/ggrastr')

library(ggrastr)

# Pass TRUE if you want to see progress output on some of Monocle 3's operations
DelayedArray:::set_verbose_block_processing(TRUE)

# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size=1000e6)

# set seed 
seed<-11979
set.seed(seed)

####################
# Helper functions
#####################

standardize <- function(z) {
  rowmed <- apply(z, 1, median)
  rowmad <- apply(z, 1, mad)  # median absolute deviation
  rv <- sweep(z, 1, rowmed,"-")  #subtracting median expression
  rv <- sweep(rv, 1, rowmad, "/")  # dividing by median absolute deviation
  return(rv)
}

lookupGeneId<-function(eset,gene_names){
  res <- rownames(fData(eset))[fData(eset)$gene_short_name %in% gene_names]
  res <- c(res,rownames(fData(eset))[rownames(fData(eset)) %in% gene_names])
  res <- unique(res)
  res
}

lookupGeneName<-function(eset,gene_id){
  res <- fData(eset[gene_id,])$gene_short_name
  res <- unique(res)
  res
}

meltCDS<-function(cds,geneset,logMode=F){
  sub<-cds[lookupGeneId(cds,geneset),]
  sub.expr<-as.matrix(exprs(sub))
  if(logMode){
    sub.expr<-log10(sub.expr+1)
  }
  sub.expr.melt<-melt(sub.expr)
  colnames(sub.expr.melt)<-c("gene_id","cell_id","value")
  res<-merge(sub.expr.melt,pData(sub),by.x="cell_id",by.y=0)
  res<-merge(res,fData(sub),by.x="gene_id",by.y=0)
  #print(head(res))
  res
}

plot_genes_barplot<-function (cds_subset, grouping = "CellType", min_expr = NULL, cell_size = 0.75, 
                              nrow = NULL, ncol = 1, panel_order = NULL, color_by = NULL, 
                              plot_trend = TRUE, label_by_short_name = TRUE, relative_expr = TRUE) 
{
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
                                                 "negbinomial.size")) {
    integer_expression <- TRUE
  }
  else {
    integer_expression <- FALSE
    relative_expr <- TRUE
  }
  if (integer_expression) {
    cds_exprs <- exprs(cds_subset)
    if (relative_expr) {
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    }
    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  }
  else {
    cds_exprs <- exprs(cds_subset)
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  }
  if (is.null(min_expr)) {
    min_expr <- cds_subset@lowerDetectionLimit
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- cds_exprs$gene_short_name
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  if (is.null(panel_order) == FALSE) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label, 
                                      levels = panel_order)
  }
  q <- ggplot(aes_string(x = grouping, y = "expression"), data = cds_exprs)
  if (plot_trend == TRUE) {
    q <- q + stat_summary(aes_string(color = color_by, fill=color_by), fun.data = "mean_cl_boot", 
                          size = 0.35, geom="bar")
    q <- q + stat_summary(aes_string(x = grouping, y = "expression", 
                                     color = color_by, group = color_by), fun.data = "mean_cl_boot", 
                          size = 0.35, geom = "errorbar")
  }
  q <- q + facet_wrap(~feature_label, nrow = nrow, 
                      ncol = ncol, scales = "free_y")
  #if (min_expr < 1) {
  #    q <- q + expand_limits(y = c(min_expr, 1))
  #}
  q <- q + ylab("Expression") + xlab(grouping)
  q <- q + monocle:::monocle_theme_opts() + theme(axis.text.x= element_text(angle = 90, hjust = 1, vjust = 0.5))
  q
}

gene_summary_heatmap<-function(cds,markers, grouping = "CellType",scale=TRUE, logMode=FALSE){
  marker.subset<-meltCDS(cds,markers,logMode = logMode)
  marker.celltype.summary<- marker.subset %>%
    group_by_(grouping,"gene_short_name") %>%
    summarise(mean_cpct=mean(value,na.rm=T))
  #print(head(marker.celltype.summary))
  
  marker.celltype.matrix<-dcast(marker.celltype.summary,as.formula(paste(grouping,"~gene_short_name")),fun.aggregate=sum)
  #print(head(marker.celltype.matrix))
  
  rownames(marker.celltype.matrix)<-marker.celltype.matrix[,1]
  marker.celltype.matrix<-marker.celltype.matrix[,-1]
  
  if(scale==TRUE){
    marker.celltype.matrix<-scale(marker.celltype.matrix)
    #print(head(marker.celltype.matrix.scaled))
  }
  #print("Ok to here!")
  marker.col.dendro<-as.dendrogram(hclust(dist(t(marker.celltype.matrix))))
  marker.row.dendro<-as.dendrogram(hclust(dist(marker.celltype.matrix)))
  
  marker.celltype.heatdata<-melt(as.matrix(marker.celltype.matrix))
  colnames(marker.celltype.heatdata)<-c(grouping,"gene_short_name","znorm")
  
  #print(head(marker.celltype.heatdata))
  #print(rownames(marker.celltype.matrix))
  
  marker.celltype.heatdata[[grouping]]<-factor(marker.celltype.heatdata[[grouping]],levels=rownames(marker.celltype.matrix)[order.dendrogram(marker.row.dendro)])
  marker.celltype.heatdata$gene_short_name<-factor(marker.celltype.heatdata$gene_short_name,levels=colnames(marker.celltype.matrix)[order.dendrogram(marker.col.dendro)])
  
  p<-ggplot(marker.celltype.heatdata) +
    geom_tile(aes_string(x='gene_short_name',y=grouping,fill='znorm')) +
    scale_fill_viridis() +
    coord_equal() + 
    theme(axis.text.x = element_text(angle=-90,hjust=0)) +
    monocle:::monocle_theme_opts()
  
  p
  
}


gene_logfc_jitter<-function(cds,markers, group = "CellType",scale=TRUE,limit_cutoff=10){
  marker.subset<-meltCDS(cds,markers)
  myCols<-unique(marker.subset)
  marker.celltype.summary<- marker.subset %>%
    group_by_(group,"Genotype","gene_short_name") %>%
    summarise(mean_cpct=mean(value,na.rm=T)) %>%
    spread(Genotype,mean_cpct) %>%
    mutate(logfc=log2(`C9-BAC`/WT))
  #print(head(marker.celltype.summary))
  myLim<-max(abs(marker.celltype.summary$logfc),na.rm=TRUE)
  p<-ggplot(marker.celltype.summary) +
    geom_jitter(aes_string(x=group,y='logfc'),color="grey80",width=0.25,data=subset(marker.celltype.summary,logfc<1 & logfc>(-1))) +
    geom_jitter(aes_string(x=group,y='logfc',color="CellType"),width=0.25,data=subset(marker.celltype.summary,logfc>=1 | logfc<=(-1))) +
    geom_hline(aes(yintercept=0),linetype="dashed") +
    theme(axis.text.x = element_text(angle=-90,hjust=0)) +
    scale_y_continuous(limits = c(-min(myLim,5), min(myLim,5))) +
    guides(color=FALSE) +
    ylab("Log2 FC (C9-BAC/WT)") + 
    monocle:::monocle_theme_opts()
  p
  
}


gene_logfc_easy_jitter<-function(cds,markers,limit_cutoff=10,colLabels=c('gene_id','gene_short_name')){
  marker.subset<-fData(cds)[lookupGeneId(cds,markers),colLabels]
  #print(head(marker.subset))
  marker.celltype.summary<-melt(as.matrix(marker.subset))
  colnames(marker.celltype.summary)<-c("gene_id","CellType","log2fc")
  #print(marker.celltype.summary)
  
  myLim<-max(abs(marker.celltype.summary$log2fc),na.rm=TRUE)
  #print(myLim)
  p<-ggplot(marker.celltype.summary) +
    geom_jitter(aes(x=CellType,y=log2fc),color="grey80",width=0.25,data=subset(marker.celltype.summary,log2fc<1 & log2fc>(-1))) +
    geom_jitter(aes(x=CellType,y=log2fc,color=CellType),width=0.25,data=subset(marker.celltype.summary,log2fc>=1 | log2fc<=(-1))) +
    geom_hline(aes(yintercept=1),linetype="dashed",color="grey60") +
    geom_hline(aes(yintercept=-1),linetype="dashed",color="grey60") +
    geom_hline(aes(yintercept=0)) +
    theme(axis.text.x = element_text(angle=-90,hjust=0)) +
    scale_y_continuous(limits = c(-min(myLim,limit_cutoff), min(myLim,limit_cutoff))) +
    guides(color=FALSE) +
    ylab("Log2 FC (C9-BAC/WT)") + 
    monocle:::monocle_theme_opts()
  p
  
}


gene_logfc_easy_heatmap<-function(cds,markers,limit_cutoff=2,colLabels=colLabels,celltypes=CellTypes){
  marker.subset<-fData(cds)[lookupGeneId(cds,markers),c('gene_short_name',colLabels)]
  rownames(marker.subset)<-make.unique(marker.subset$gene_short_name)
  marker.subset<-marker.subset[,-1]
  colnames(marker.subset)<-celltypes
  marker.subset[is.na(marker.subset)]<-0
  
  myLim<-min(limit_cutoff,max(abs(marker.subset),na.rm=TRUE))
  marker.subset[marker.subset>myLim]<-myLim
  marker.subset[marker.subset<(-myLim)]<-(-myLim)
  #print(head(marker.subset))
  
  marker.celltype.summary<-melt(as.matrix(marker.subset))
  colnames(marker.celltype.summary)<-c("gene_id","CellType","log2fc")
  #print(marker.celltype.summary)
  
  marker.celltype.summary$log2fc[marker.celltype.summary$log2fc>myLim]<-myLim
  marker.celltype.summary$log2fc[marker.celltype.summary$log2fc<(-myLim)]<-(-myLim)
  marker.celltype.summary$log2fc[is.na(marker.celltype.summary$log2fc)]<-0
  
  row.dendro<-order.dendrogram(as.dendrogram(hclust(dist(marker.subset))))
  col.dendro<-order.dendrogram(as.dendrogram(hclust(dist(t(marker.subset)))))
  #print(rownames(marker.subset)[row.dendro])
  
  marker.celltype.summary$gene_id<-factor(marker.celltype.summary$gene_id,levels=rownames(marker.subset)[row.dendro])
  marker.celltype.summary$CellType<-factor(marker.celltype.summary$CellType,levels=colnames(marker.subset)[col.dendro])
  #print(myLim)
  p<-ggplot(marker.celltype.summary) +
    geom_tile(aes(x=CellType,y=gene_id,fill=log2fc)) +
    theme(axis.text.x = element_text(angle=-90,hjust=0)) +
    scale_fill_gradient2(low="steelblue",mid="white",high="red") +
    coord_equal(0.3) +
    monocle:::monocle_theme_opts()
  p
  
}


