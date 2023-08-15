#' Annotate complex event blocks
#'
#' This function annotates complex event blocks by identifying genes that overlap with the specified genomic regions. It utilizes the GenomicRanges package to perform overlap detection and retrieves gene symbols for the overlapping genes.
#'
#' @param block A data frame or tibble containing complex event blocks. Each row should represent a block with "chr" (chromosome) and "start" (start position) columns.
#' @param window_size The window size used for creating genomic regions around blocks (default is 10000).
#' @param USCS_gene The name of the UCSC transcript database object (default is "TxDb.Hsapiens.UCSC.hg38.knownGene").
#' @return A character vector containing unique gene symbols associated with overlapping blocks.
#' @importFrom tidyr separate
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicFeatures genes
#' @importFrom org.Hs.eg.db select
#' @importFrom AnnotationDbi select
#' @importFrom base gsub
#' @importFrom base as.numeric
#' @export

complex_block_annotation <- function(block,window_size=10000,USCS_gene="TxDb.Hsapiens.UCSC.hg38.knownGene") {
  library(tidyr)
  library(GenomicRanges)
  library("GenomicFeatures")
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  txdb.gene <-genes(txdb, columns=c("TXCHROM","TXSTART", "TXEND","GENEID","TXNAME"))
  blocks <- as.data.frame(block)
  colnames(blocks) <- "blocks"
  blocks_tmp <- separate(blocks, blocks, c("chr", "start"))
  blocks_tmp$chr <- gsub("C", "chr", blocks_tmp$chr)
  blocks_tmp$start <- as.numeric(blocks_tmp$start)
  blocks_tmp$start <- blocks_tmp$start * window_size
  blocks_tmp$end <- blocks_tmp$start + window_size
  blocks_tmp_GR <- makeGRangesFromDataFrame(blocks_tmp,
                                            keep.extra.columns = TRUE,
                                            ignore.strand = FALSE)
  hits <- findOverlaps(txdb.gene, blocks_tmp_GR)
  txdb.gene_hits <- txdb.gene[queryHits(hits), ]
  ENTREZID_gene <- AnnotationDbi::select(
                   org.Hs.eg.db,
                   keys = as.character(txdb.gene_hits$GENEID),
                   columns = c("SYMBOL")
                   )
  block_gene <- unique(ENTREZID_gene$SYMBOL)
  return (block_gene)
}
