library(Pirat)

#' Remove peptides or ions with at least t MVs
#'
#' @param l_pep_rna 
#' @param t 
#' @param percentage 
#'
#' @return
#' @export
#'
#' @examples
remove_NA_pep_capizzi <- function(l_pep_rna, t, percentage=F) {
  nsample = nrow(l_pep_rna$peptides_ab)
  if (percentage) {
    t = ceiling(nsample*t)
  }
  mv_count_pep = colSums(is.na(l_pep_rna$peptides_ab))
  i_pep_rm = which(mv_count_pep >= t)
  if (length(i_pep_rm) != 0) {
    new_adj = l_pep_rna$adj[-i_pep_rm,]
    i_rna_rm = which(colSums(new_adj) == 0)
    new_pep = l_pep_rna$peptides_ab[,-i_pep_rm]
    if (length(i_rna_rm) != 0) {
      new_adj = new_adj[,-i_rna_rm]
      return(list("peptides_ab" = new_pep, "adj" = new_adj))
    }
    else {
      return(list("peptides_ab" = new_pep,
                  "adj" = new_adj))
    }
  }
  else {
    return(l_pep_rna)
  }
}

#' Data loader for Capizzi2022
#'
#' @param path2data 
#' @param n.na.max 
#' @param rm.nested.pg 
#'
#' @return
#' @export
#'
#' @examples
get_capizzi_data = function(path2data, n.na.max = NULL, rm.nested.pg = F) {
  peptide_data = read.csv(path2data, sep = ";")
  peptide_data = peptide_data[peptide_data$Status == "Validated automatically", ]
  allproteins = as.character(peptide_data$Protein.Sets)
  alluniq_pep = paste(peptide_data$Peptide.Sequence, peptide_data$PTMs)
  
  samplesname = c("Abundance.WT.5", "Abundance.WT.4",
                  "Abundance.WT.3", "Abundance.WT.2",
                  "Abundance.WT.1", "Abundance.CAG.5",
                  "Abundance.CAG.4", "Abundance.CAG.3",
                  "Abundance.CAG.2", "Abundance.CAG.1"
  )
  
  mask_common_prots = grepl(", ", allproteins)
  common_prots = as.character(allproteins[mask_common_prots])
  all_common_prots = unlist(sapply(common_prots, strsplit, ", ")) # Get prots from shared peptides
  alluniq_prots = unique(c(unique(all_common_prots), unique(allproteins[!mask_common_prots])))

  pep.data = matrix(t(peptide_data[, samplesname]), nrow = length(samplesname))
  colnames(pep.data) = alluniq_pep

  adj = matrix(F, ncol = length(alluniq_prots), nrow = length(alluniq_pep))
  colnames(adj) = alluniq_prots
  rownames(adj) = alluniq_pep

  for (i in 1:length(alluniq_prots)) {
    adj[, i] = grepl(paste(alluniq_prots[i], "($|, )", sep = ""), allproteins)
  }
  print(table(rowSums(adj)))

  # Remove contaminants
  mask_con_prot = grepl("#C#", alluniq_prots)
  adj = adj[, !mask_con_prot]
  mask_con_pep = rowSums(adj) == 0
  pep.data = pep.data[, !mask_con_pep]
  adj = adj[!mask_con_pep, ]

  pep.data[pep.data == 0] = NA

  processed.data = list(peptides_ab = t(vsn::justvsn(t(pep.data))), adj = adj)

  if (!is.null(n.na.max)) {
    processed.data = remove_NA_pep_capizzi(processed.data, n.na.max)
  }

  if (rm.nested.pg) {
    idx.emb.prots = get_indexes_embedded_prots(processed.data$adj)
    processed.data = rm_pg_from_idx_merge_pg(processed.data, idx.emb.prots)
  }

  return(processed.data)
}



