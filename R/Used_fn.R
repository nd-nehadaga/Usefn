##### To make gr object from location coordinate
#' Create unique entry gr object with metacolumns
#'
#' @param genomic_coord Vector of genomic coordinates
#' @param metacolumn_name Genomic coordinate's column name
#'
#' @return Gr object with unique entries
#' @export
#'
#' @examples unique_grobject_metacol(genomic_coord = c("chr2:15000-15500","chr5:10001-10122"), metacolumn_name = "Peak_id")
#' @examples unique_grobject_metacol(genomic_coord = c("chr2:15000-15500","chr5:10001-10122","chr2:15000-15500"), metacolumn_name = "Peak_id")

unique_grobject_metacol = function (genomic_coord,metacolumn_name) {

  chr_loc = unique(genomic_coord)
  chr_loc_gr = GenomicRanges::GRanges(chr_loc)
  chr_loc_gr$meta_col = chr_loc

  names(GenomicRanges::mcols(chr_loc_gr)) = metacolumn_name
  return (chr_loc_gr)  }



#' Overlap between two genomic range objects with length of intersection
#'
#' @param gr1 Genomic Range Object1
#' @param gr2 Genomic Range Object2
#' @param gr1_col Column name for gr1
#' @param gr2_col Column name for gr2
#'
#' @return Data frame with overlap and length of intersection
#' @export
#'
#' @examples
#'

to_overlap_grange_object = function (gr1, gr2, gr1_col, gr2_col) {

  query_gr = gr1
  subject_gr = gr2

  query_col = gr1_col
  subject_col = gr2_col

  overlap = IRanges::findOverlaps(query_gr,subject_gr)

  query_row_ids = queryHits(overlap)
  subject_rowids = subjectHits(overlap)

  query_overlap_df = as.data.frame(GenomicRanges::elementMetadata(query_gr)[query_row_ids,query_col])
  subject_overlap_df= as.data.frame(GenomicRanges::elementMetadata(subject_gr)[subject_rowids, subject_col])

  ### to find the interesection of overlaps
  overlap_width = width(pintersect(query_gr[query_row_ids],subject_gr[subject_rowids]))

  Overlap_df = cbind.data.frame(query_overlap_df,subject_overlap_df)
  names(Overlap_df) = c(gr1_col, gr2_col)
  Overlap_df["width"] = overlap_width



  return(Overlap_df)

}

#' Nearest distance grane objects with the distance value
#'
#' @param gr1 Genomic Range Object1
#' @param gr2 Genomic Range Object2
#' @param gr1_col Column name for gr1
#' @param gr2_col Column name for gr2
#'
#' @return Data frame with genomic range objects with nearest distance
#' @export
#'
#' @examples


to_find_nearest_dist_grange_object = function (gr1, gr2, gr1_col, gr2_col) {

  query_gr = gr1
  subject_gr = gr2

  query_col = gr1_col
  subject_col = gr2_col

  overlap = IRanges::distanceToNearest(query_gr,subject_gr,select="arbitrary")

  query_row_ids = queryHits(overlap)
  subject_rowids = subjectHits(overlap)

  query_overlap_df = as.data.frame(GenomicRanges::elementMetadata(query_gr)[query_row_ids,query_col])
  subject_overlap_df= as.data.frame(GenomicRanges::elementMetadata(subject_gr)[subject_rowids, subject_col])
  dist_df = as.data.frame(GenomicRanges::elementMetadata(overlap))

  Overlap_df = Reduce(f = cbind.data.frame,list(query_overlap_df,subject_overlap_df,dist_df))

  return(Overlap_df)

}



#' Overlap between two sets of genomic coordinates as input
#'
#' @param input_genomic_coord1 Vector1 containing genomic coordinates
#' @param input_genomic_coord2 Vector2 containing genomic coordinates
#' @param input1_colname Column name of input vector1
#' @param input2_colname Column name of input vector2
#'
#' @return Data frame with overlap and length of intersection
#' @export
#'
#' @examples

to_get_overlap_bw_gr_from_peakids = function (input_genomic_coord1,input_genomic_coord2,input1_colname,input2_colname) {

  input_genomic_coord1 = input_genomic_coord1[!(is.na(input_genomic_coord1))]
  input_genomic_coord2 = input_genomic_coord2[!(is.na(input_genomic_coord2))]


  input1_gr = unique_grobject_metacol(genomic_coord = input_genomic_coord1,metacolumn_name = "Location")
  input2_gr = unique_grobject_metacol(genomic_coord = input_genomic_coord2,metacolumn_name= "Location")


  overlap_df = to_overlap_grange_object(gr1 = input1_gr,gr2 = input2_gr,
                                        gr1_col = "Location",gr2_col = "Location")

  names(overlap_df) = c(input1_colname,input2_colname,"width")

  return(overlap_df)

}

############ Modify strings



#' Split string based on separator
#'
#' @param x Input string
#' @param delimtter Splitting separator
#' @param position Position of output string after splitting
#'
#' @return Split string
#' @export
#'
#' @examples


modify_string = function (x,delimtter,position) {


  output_string = unlist(strsplit(as.character(x),split = delimtter))[position]

  return (output_string)

}

###### To get bed format from peak coordinate (location)

#' Convert genomic coordinates into bed format
#'
#' @param input_peak_id_vector Vector with genomic coordinates
#' @param bed_file Default parameter is FALSE
#' @param out_fn Output file name for writing bed format. Deafault is "
#'
#' @return Data frame and output file with genomic coordinates in bed format
#' @export
#'
#' @examples


to_get_bed_df_peaks_from_peak_id_vector = function(input_peak_id_vector,bed_file = FALSE, out_fn = "") {

  chr_name = sapply(X = input_peak_id_vector,FUN = modify_string,delimtter = ":",position = 1)
  chr_str_stop = sapply(X = input_peak_id_vector,FUN = modify_string,delimtter = ":",position = 2)
  chr_str = sapply(X = chr_str_stop,FUN = modify_string,delimtter = "-",position = 1,USE.NAMES = F)
  chr_stop = sapply(X = chr_str_stop,FUN = modify_string,delimtter = "-",position = 2,USE.NAMES = F)

  peak_id_df = data.frame("Chr_name" = chr_name, "start" = chr_str, "stop" = chr_stop)

  if (bed_file) {

    write.table(x = peak_id_df,file = out_fn,sep = "\t",quote = F,row.names = F,col.names = F)

  }

  return(peak_id_df)
}


################### PCA, heatmap,correlogram, ggpairs
#' Computing PCA
#'
#' @param data_mat Count matrix rld or vst transformed and normalised
#' @param ntop Top n variable features
#' @param meta_df Data frame of covariates
#'
#' @return Data frame with PCA loadings and covariates
#' @export
#'
#' @examples



compute_pca <- function(data_mat, ntop,meta_df){

  Pvars <- rowVars(data_mat)
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop,
                                                        length(Pvars)))]

  PCA <- prcomp(t(data_mat[select, ]), scale. = F, center = TRUE)
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

  pca_df = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                      PC3 = PCA$x[,3], PC4 = PCA$x[,4],PC5 = PCA$x[,5],PC6 = PCA$x[,6])

  #meta_df = as.data.frame(colData(data_mat))

  pca_meta_df = cbind(pca_df, meta_df)

  return(list(PCA_meta_df = pca_meta_df,
              perc_var = percentVar,
              selected_vars = select))

}


######## Correlogram
#' Correlogram of pca loadings and covariates
#'
#' @param count_matrix Count matrix rld or vst transformed normalised
#' @param ntop top n features
#' @param metadata Covariates data frame
#' @param covs Column names of covariates
#' @param scale logical value TRUE or FALSE
#' @param pcCount Number of dimension
#' @param outdir Output file name for correlogram
#' @param logTransform Deafault is TRUE
#'
#' @return Heatmap of correlogram between PCA loadings and covariates
#' @export
#'
#' @examples

pca.correlogram <- function(count_matrix,ntop, metadata, covs, scale, pcCount, outdir, logTransform = TRUE){


  Pvars <- rowVars(count_matrix)
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop,
                                                        length(Pvars)))]



  ## adj. R2-values ##

  r2_fun <- function(x, y) {
    summary(lm(y ~ x))$adj.r.squared
  } #function to get adj. R2 values

  PCA <- prcomp(t(count_matrix[select, ]), scale = scale)
  pcs_cv <- PCA$x[,1:pcCount]   #Top PCs from top 1000 variable genes
  eig.val <- 100*summary(PCA)$importance[2,1:pcCount]
  labels <- paste0(colnames(pcs_cv), "_",  "(", round(eig.val, 2),  " %)")

  ###### Metadata

  covariates <- data.frame(metadata[,covs],
                           row.names=colnames(count_matrix))    #Clinical covariates

  r2_cv <- matrix(NA, nrow = ncol(covariates), ncol = ncol(pcs_cv),
                  dimnames = list(colnames(covariates), colnames(pcs_cv)))

  for (cov in colnames(covariates)) {
    for (pc in colnames(pcs_cv)) {
      r2_cv[cov, pc] <- r2_fun(covariates[,cov], pcs_cv[,pc])
    }
  }      #get adj. R2 value for covariates ~ PC 1-6

  r2CvPlot <- pheatmap::pheatmap(r2_cv, cluster_cols = F, display_numbers = T,labels_col = labels,
                       color = colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = "Reds")))(1000))

  xlabCur = "Mean log counts"

  if (logTransform == TRUE) {
    densPlot <- function() {
      dens <- apply(log2(count_matrix+1), 2, function(x) density(x, na.rm = TRUE))
      plot(NA,xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")))
      mapply(lines, dens, col=1:length(dens))
    }
  }
  else {
    densPlot <- function() {
      dens <- apply(count_matrix, 2, density)

      plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")))
      mapply(lines, dens, col=1:length(dens))

    }
  }

  pdf(outdir, height=10, width=10)

  grid::grid.draw(r2CvPlot$gtable)
  densPlot()

  dev.off()

}


#' Clustered heatmap based on variable features
#'
#' @param data_mat Count matrix (rld or vst transformed, normalised)
#' @param ntop Top n variable features
#' @param col_select Vector of annotation column names
#' @param meta_df Row annotation data frame
#' @param outfn Output filename
#' @param rownorm Z-score across rows. Logical value TRUE/FALSE
#'
#' @return Heatmap based on n variable features
#' @export
#'
#' @examples

draw_heatmap_top_peaks = function (data_mat,ntop,col_select,meta_df,outfn,rownorm) {

  Pvars <- rowVars(data_mat)
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop,
                                                        length(Pvars)))]

  ##### DataFrame
  input_df = as.data.frame(data_mat)
  input_df$rowvar = Pvars
  input_df = input_df[order(input_df$rowvar,decreasing = TRUE),]

  ####### Top peak mat
  data_mat_sub = data_mat[select, ]


  ######### Normalises across rowa
  cal_z_score <- function(x){
    (x - mean(x,na.rm = TRUE)) / sd(x,na.rm = TRUE)
  }

  data_mat_norm <- t(apply(data_mat_sub, 1, cal_z_score))


  ############# To get input mat and if row normalised or not
  if (rownorm == TRUE) {
    input_mat = data_mat_norm

  }else{
    input_mat = data_mat_sub
  }


  ##### For annotation
  #meta_df = as.data.frame(colData(data_mat))
  annotate_col_df = meta_df[,col_select]
  rownames(annotate_col_df) = rownames(meta_df)

  annotate_row_df = data.frame("Peak_var" = input_df$rowvar[1:ntop])
  rownames(annotate_row_df) = rownames(input_df)[1:ntop]

  ########Check
  check1 = all(rownames(annotate_row_df) == rownames(input_mat))
  check2 = all(colnames(input_mat) == rownames(annotate_col_df))


  if (check1 & check2) {

    p1 = pheatmap::pheatmap(input_mat,annotation_col = annotate_col_df,cluster_rows = T,cluster_cols = T,show_rownames = F,annotation_row = annotate_row_df)
    p2 = pheatmap::pheatmap(input_mat,annotation_col = annotate_col_df,cluster_rows = F,cluster_cols = F,show_rownames = F,annotation_row = annotate_row_df)

    pdf(file = outfn,width = 10,height = 10,onefile = TRUE,useDingbats = FALSE)

    grid::grid.newpage()
    grid::grid.draw(p1$gtable)

    grid::grid.newpage()
    grid::grid.draw(p2$gtable)

    dev.off()

  }

}



##### Diagniostic plots related to differential expression (MA, volcano plot)

#' MA Plot labelled with significant genes
#'
#' @param res DESeq2 result table of differential analysis
#' @param thresh Adjusted p-value threshold. Default is 0.05
#' @param labelsig Label genes below threshold. Logical value TRUE or FALSE. Default is TRUE
#' @param textcx Text size. Default is 1
#' @param ...Other. Parameters for base R plot function
#'
#' @return MA plot for differential analysis
#' @export
#'
#' @examples


maplot <- function (res, thresh=0.05, labelsig=FALSE, textcx=1, ...) {

  res = as.data.frame(res)
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))

  if (labelsig) {
    require(calibrate)
    Gene = rownames(subset(res, padj<thresh))
    base::with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}



## Volcano plot with "significant" genes labeled

#' Volcano plot with "significant" genes labeled
#'
#' @param res DESeq2 result table of differential analysis
#' @param lfcthresh Log2 Fold change threshold. Default is 2
#' @param sigthresh Significant p-value threshold. Default is 0.05
#' @param main plot title
#' @param legendpos Position of legend . The value can be c()
#' @param labelsig Label genes below sigthresh and above lfcthresh. Logical value TRUE or FALSE. Default is TRUE
#' @param textcx text size. Default is 1
#' @param ...Parameters other parameters for base R plot function
#'
#' @return Volcano plot with labelled genes
#'
#' @export
#'
#' @examples

volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topright", labelsig=FALSE, textcx=1, ...) {

  res = as.data.frame(res)
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    Gene = rownames(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh))
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}




##### Mapping gene ids
#' Mapping ENSEMBL to gene symbols
#'
#' @param input_ensembl_vector Vector of ensembl ids to be mapped
#'
#' @return Data frame with gene symbols, entrezid
#'
#' @export
#'
#' @examples

to_add_gene_symbol_ENSEMBL = function (input_ensembl_vector) {

  ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

  gene_mapping_df = biomaRt::getBM(mart = ensembl,attributes = c("entrezgene_id","ensembl_gene_id","hgnc_symbol"),
                          filters = "ensembl_gene_id",values = input_ensembl_vector)
  return (gene_mapping_df)
}



###### Enrichment analysis
#' Title
#'
#' @param expressedTarget Forground genes (ENTREZ IDs)
#' @param expressedGenes Background/universe genes (ENTREZ IDs)
#' @param pval_thresh Significant p-val threshold. Default is 0.05
#' @param ont Default is biological process (BP). Value can be any/all GO term "BP", "MF" or "CC
#'
#' @return List of GO and kegg enrichment
#'
#' @export
#'
#' @examples


to_get_GO_kegg_enrichment = function (expressedTarget,expressedGenes,pval_thresh = 0.05,ont = "BP") {


  GO <- clusterProfiler::enrichGO(gene = expressedTarget,
                 OrgDb='org.Hs.eg.db',
                 pvalueCutoff= 1,
                 ont = ont,
                 keyType= "ENTREZID",
                 universe=expressedGenes)

  GO_df= GO@result
  GO_df = GO_df[order(GO_df$p.adjust),]

  kegg <- clusterProfiler::enrichKEGG(gene = expressedTarget,
                   organism  = "human",
                   pvalueCutoff = 0.05,universe = expressedGenes)

  kegg_df =  kegg@result
  kegg_df =  kegg_df[order(kegg_df$p.adjust),]


  return (list(GO_df,kegg_df))

}


#' Fisher exact test
#'
#' @param enrichment_category_groups Functional enrichment group name vector (colnames)
#' @param input_category_groups_to_be_tested Input category group name vector (rownames)
#' @param x overlap of hit list with category of interest
#' @param n total no of genes (universe)
#' @param m  number of genes in a category of interest
#' @param k length of hit list
#' @param ... other parameters to be passed in fisher.test function
#'
#' @return List of contingency table and result of fisher test
#'
#' @export
#'
#' @examples

to_do_fisher_exact_test = function (enrichment_category_groups,input_category_groups_to_be_tested,x,n,m,k,...) {

  #k = length of hit list;
  #x = overlap of hit list with category of interest;
  #m = number of genes in a category of interest;
  #n = total no of genes (universe)

  data_vector = c(x,k-x,m-x,((n-m)-(k-x)) )

  contingency.table = matrix(data_vector,nrow = 2,ncol = 2,byrow = TRUE)
  colnames(contingency.table) <- enrichment_category_groups
  rownames(contingency.table) <- input_category_groups_to_be_tested


  f.test = fisher.test(contingency.table,...)

  return (list(f.test,contingency.table))

}



