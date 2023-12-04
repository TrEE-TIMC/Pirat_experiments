source("../method/utils.R")

#' Remove peptides or ions with at least t MVs on datasets with Ground truth DA peptides
#'
#' @param l_pep_rna 
#' @param t 
#' @param percentage 
#'
#' @return
#' @export
#'
#' @examples
remove_NA_pep_maxquant <- function(l_pep_rna, t, percentage=F) {
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
    new_mask_pep_diff = l_pep_rna$mask_pep_diff[-i_pep_rm]
    if (length(i_rna_rm) != 0) {
      new_mask_prot_diff = l_pep_rna$mask_prot_diff[-i_rna_rm]
      new_adj = new_adj[,-i_rna_rm]
      return(list("peptides_ab" = new_pep, "adj" = new_adj,
                  mask_prot_diff = new_mask_prot_diff,
                  mask_pep_diff = new_mask_pep_diff))
    }
    else {
      return(list("peptides_ab" = new_pep,
                  "adj" = new_adj, mask_prot_diff = l_pep_rna$mask_prot_diff,
                  mask_pep_diff = new_mask_pep_diff))
    }
    
  }
  else {
    return(l_pep_rna)
  }
}

#' Title
#'
#' @param peptide_path 
#' @param suppmat2_path 
#'
#' @return
#' @export
#'
#' @examples
get.maxquant.data = function(peptide_path, suppmat2_path, n.na.max = NULL, rm.nested.pg = F) {
  peptide_data = read.csv(peptide_path, sep="\t")
  allrazorprot = peptide_data$Leading.razor.protein
  maskcon = grepl("CON__", allrazorprot) | grepl("REV__", allrazorprot) # Remove contaminants and REV ??
  peptide_data = peptide_data[!maskcon,]
  samplesname = c("Intensity.UPS1_01", "Intensity.UPS1_02", "Intensity.UPS1_03",
                  "Intensity.UPS1_04", "Intensity.UPS2_01", "Intensity.UPS2_02",
                  "Intensity.UPS2_03", "Intensity.UPS2_04")
  
  #Extract peptides ab
  peptides_ab = t(peptide_data[,samplesname])
  rownames(peptides_ab) = samplesname
  colnames(peptides_ab) = peptide_data$Sequence
  peptides_ab[peptides_ab == 0] = NA
  
  allproteins = peptide_data$Proteins
  mask_common_prots = grepl(";", allproteins)
  common_prots = allproteins[mask_common_prots]
  all_common_prots = unlist(sapply(common_prots, strsplit, ";")) # Get prots from shared peptides
  uniq_prots = unique(c(all_common_prots, allproteins[!mask_common_prots]))
  uniq_prots = uniq_prots[!grepl("CON__", uniq_prots)]
  adj = matrix(0, nrow = ncol(peptides_ab), ncol = length(uniq_prots))
  colnames(adj) = uniq_prots
  rownames(adj) = peptide_data$Sequence
  for (i in 1:length(uniq_prots)) {
    adj[, i] = as.numeric(grepl(uniq_prots[i], allproteins))
  }
  protein_data = read.csv(suppmat2_path, sep=';', dec=",")
  upsproteins = protein_data[grepl("ups", protein_data$Majority.protein.IDs), ]
  upsproteins_diff = upsproteins[upsproteins$UPS2.1.ratio != 1, ]
  mask_prot_diff = colSums(sapply(uniq_prots, grepl, upsproteins_diff$Majority.protein.IDs)) >= 1
  mask_pep_diff = rowSums(adj[, mask_prot_diff]) >= 1

  processed.data = list(peptides_ab = log(peptides_ab, 2), adj = adj, 
                        mask_prot_diff = mask_prot_diff, mask_pep_diff = mask_pep_diff)
  
  if (!is.null(n.na.max)) {
    processed.data = remove_NA_pep_maxquant(processed.data, n.na.max)
  }
  
  if (rm.nested.pg) {
    idx.emb.prots = get_indexes_embedded_prots(processed.data$adj)
    processed.data = rm_pg_from_idx_merge_pg(processed.data, idx.emb.prots)
  }
  
  return(processed.data)
}

# # Remove shared peptides
# remove_shared_pep_maxquant <- function(l_pep_rna) {
#   nsample = nrow(l_pep_rna$peptides_ab)
#   i_pep_rm = which(rowSums(l_pep_rna$adj) >= 2)
#   if (length(i_pep_rm) != 0) {
#     new_adj = l_pep_rna$adj[-i_pep_rm,]
#     i_rna_rm = which(colSums(new_adj) == 0)
#     new_pep = l_pep_rna$peptides_ab[,-i_pep_rm]
#     new_mask_pep_diff = l_pep_rna$mask_pep_diff[-i_pep_rm]
#     if (length(i_rna_rm) != 0) {
#       new_mask_prot_diff = l_pep_rna$mask_prot_diff[-i_rna_rm]
#       new_adj = new_adj[,-i_rna_rm]
#       return(list("peptides_ab" = new_pep, "adj" = new_adj,
#                   mask_prot_diff = new_mask_prot_diff,
#                   mask_pep_diff = new_mask_pep_diff))
#     }
#     else {
#       return(list("peptides_ab" = new_pep,
#                   "adj" = new_adj, mask_prot_diff = l_pep_rna$mask_prot_diff,
#                   mask_pep_diff = new_mask_pep_diff))
#     }
#   }
#   else {
#     return(l_pep_rna)
#   }
# }