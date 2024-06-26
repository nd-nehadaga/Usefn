% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Used_fn.R
\name{to_get_bed_df_peaks_from_peak_id_vector}
\alias{to_get_bed_df_peaks_from_peak_id_vector}
\title{Convert genomic coordinates into bed format}
\usage{
to_get_bed_df_peaks_from_peak_id_vector(
  input_peak_id_vector,
  bed_file = FALSE,
  out_fn = ""
)
}
\arguments{
\item{input_peak_id_vector}{Vector with genomic coordinates}

\item{bed_file}{Default parameter is FALSE}

\item{out_fn}{Output file name for writing bed format. Deafault is "}
}
\value{
Data frame and output file with genomic coordinates in bed format
}
\description{
Convert genomic coordinates into bed format
}
