#' Helper function to convert_probe_to_gene
#'
#' @param probe_ids list of probe ids
#' @param hm450 Large GRanges object
#'
#' @return vector of mapped genes
#' @export
#'
#' @examples
#' library(FDb.InfiniumMethylation.hg19)
#' hm450 <- FDb.InfiniumMethylation.hg19::get450k()
#' get_genes(c("cg10861017", "cg10861599", "cg10861751"), hm450)
get_genes <- function(probe_ids, hm450) {
  probes <- hm450[probe_ids]
  nearest_genes <- FDb.InfiniumMethylation.hg19::getNearestGene(probes)[4]
  # we want to return just the nearestGeneSymbols
  return(nearest_genes$nearestGeneSymbol)
}

#' Convert probe ids to nearest genes
#'
#' Uses FDb.InfiniumMethylation.hg19 Bioconductor package to map probe ids to nearest genes in chunks
#' @param probe_ids vector
#' @param icgc_donor_id vector
#' @param hm450 Large GRanges object
#' @param chunk_size size of chunks when converting probe ids to genes
#'
#' @return vector of mapped genes
#' @export
#'
#' @examples
#' library(FDb.InfiniumMethylation.hg19)
#' hm450 <- FDb.InfiniumMethylation.hg19::get450k()
#' convert_probe_to_gene(
#'   c("cg10861017", "cg10861599", "cg10861751", "cg10862535", "cg10862587"),
#'   "DO110234",
#'   hm450
#' )
convert_probe_to_gene <- function(probe_ids, icgc_donor_id, hm450=NULL, chunk_size = 10000) {

  if(is.null(hm450)){
    hm450 <- FDb.InfiniumMethylation.hg19::get450k()
  }

  probe_id_chunks <- split(probe_ids, ceiling(seq_along(probe_ids) / chunk_size))
  nearest_genes <- c()
  cat("processing donor: ", icgc_donor_id, "\n")

  i <- 1
  print("======= processing chunks ========")
  for (probe_id_chunk in probe_id_chunks) {
    print(i)
    i <- i + 1
    res <- tryCatch(get_genes(probe_id_chunk, hm450), error = function(x) {
      cat("ERROR with donor: ", icgc_donor_id)
      return(rep("ERROR", length(probe_id_chunk)))
    })
    nearest_genes <- append(nearest_genes, res)
  }
  return(nearest_genes)
}
