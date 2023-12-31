% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_junctions.r
\name{extract_junctions}
\alias{extract_junctions}
\title{Extract contig junctions}
\usage{
extract_junctions(one_complex_event, window_size = 10000, ...)
}
\arguments{
\item{one_complex_event}{A dataframe representing a single complex event with relevant columns including "mapID", "chr", "start_q", "end_q", "start_r", "end_r", and "strand".}

\item{window_size}{The window size used for creating genomic regions around contig junctions (default is 10000).}

\item{...}{Additional arguments (currently not used).}
}
\value{
A dataframe containing extracted contig junctions with combined blocks and gaps.
}
\description{
This function takes a complex event dataframe and extracts contig junctions based on specified window size and gap threshold. It identifies contig junctions within each complex event, considering strand information and block rankings.
}
