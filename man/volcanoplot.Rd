% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Used_fn.R
\name{volcanoplot}
\alias{volcanoplot}
\title{Volcano plot with "significant" genes labeled}
\usage{
volcanoplot(
  res,
  lfcthresh = 2,
  sigthresh = 0.05,
  main = "Volcano Plot",
  legendpos = "topright",
  labelsig = FALSE,
  textcx = 1,
  ...
)
}
\arguments{
\item{res}{DESeq2 result table of differential analysis}

\item{lfcthresh}{Log2 Fold change threshold. Default is 2}

\item{sigthresh}{Significant p-value threshold. Default is 0.05}

\item{main}{plot title}

\item{legendpos}{Position of legend . The value can be c()}

\item{labelsig}{Label genes below sigthresh and above lfcthresh. Logical value TRUE or FALSE. Default is TRUE}

\item{textcx}{text size. Default is 1}

\item{...Parameters}{other parameters for base R plot function}
}
\value{
Volcano plot with labelled genes
}
\description{
Volcano plot with "significant" genes labeled
}
