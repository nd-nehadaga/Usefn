#' Data to be used as examples for the function usage
#'
#' Contains epigenetic profiling H3K27ac ChIP-seq count data together with meta file for two tissues PBMC (blood) and synovial fluid (SF) in JIA patients
#' The data is publicly available
#'
#' @format Matrix with count data which is rld transformed having 26730 peaks (rows) and 8 samples (column)

#'
#' @source The raw data is present here https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71596
#'
#' @examples
#' data(count_matrix)
"count_matrix"

#' Contains metadata associated with epigenetic profiling H3K27ac ChIP-seq count data for two tissues PBMC (blood) and synovial fluid (SF) in JIA patients
#' @format Data frame with meta data which associated with count data
#' @source The raw data is present here https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71596
#' @examples
#' data(meta_df)
"meta_df"

#' Contains differential gene expression data between JIA patients (SF) and healthy controls (PBMC) . It is result table from DESeq2 analysis.
#' @format Data frame of DESeq2 differential analysis result object
#' @source The raw data is present here https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71595
#' @examples
#' data(res_df)
"res_df"
