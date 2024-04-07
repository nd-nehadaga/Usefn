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
