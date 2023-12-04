library(Pirat)

put_imputed_in_or = function(imputed.ds, or.ds) {
  imputed.names = colnames(imputed.ds)
  or.names = colnames(or.ds)
  is.imputed = or.names %in% imputed.names
  or.ds[, is.imputed] = imputed.ds
  return(or.ds)
}

#' @title Remove peptides by index
#' 
#' @description Remove peptide by index, then deletes empty PGs
#'
#' @param l_pep_rna List representing the dataset
#' @param pepidxs Vector of indices
#'
#' @export
#' @return List representing the dataset
#'
remove_pep_from_idx <- function(l_pep_rna,
                                pepidxs) {
  i_pep_rm = pepidxs
  if (length(i_pep_rm) != 0) {
    l_pep_rna$adj = l_pep_rna$adj[-i_pep_rm,]
    l_pep_rna$peptides_ab = l_pep_rna$peptides_ab[,-i_pep_rm]
    if (!is.null(l_pep_rna$mask_pep_diff)) {
      l_pep_rna$mask_pep_diff = l_pep_rna$mask_pep_diff[-i_pep_rm]
    }
    if (!is.null(l_pep_rna$charges)) {
      l_pep_rna$charges = l_pep_rna$charges[-i_pep_rm]
    }
    if (!is.null(l_pep_rna$modifs)) {
      l_pep_rna$modifs = l_pep_rna$modifs[-i_pep_rm]
    }
    i_pg_rm = which(colSums(l_pep_rna$adj) == 0)
    if (length(i_pg_rm) != 0) {
      l_pep_rna = rm_pg_from_idx_merge_pg(l_pep_rna, i_pg_rm)
    }
  }
  return(l_pep_rna)
}