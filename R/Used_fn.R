##### To make gr object from location coordinate
#' Make unique entry gr object with metacolumns
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



