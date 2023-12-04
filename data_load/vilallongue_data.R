source("../data_loader/capizzi_data.R")

#' Data loader for Vilallongue2022
#'
#' @param path2folder 
#' @param pep_ion 
#' @param target 
#' @param n.na.max 
#' @param rm.nested.pg 
#'
#' @return
#' @export
#'
#' @examples
get_vilallongue_data = function(path2folder, pep_ion = "ion", target="SC",
                                n.na.max = NULL, rm.nested.pg = F, norm.method = NULL) {
  
  path2data = file.path(path2folder, paste0(target, "_", pep_ion, ".csv"))
  peptide_data = read.csv(path2data, sep = ";")
  if (pep_ion == "peptide") {
    peptide_data = peptide_data[peptide_data$Status == "Validated automatically", ]
  }
  allproteins = as.character(peptide_data$Protein.Sets)
  
  alluniq_pep = paste(peptide_data$Peptide.Sequence, peptide_data$PTMs, sep="_")
  if (pep_ion == "ion") {
    alluniq_pep = paste(alluniq_pep, peptide_data$Charge, sep = "_")
  }
  
  if (target == "dLGN") {
    samplesname = c("Abundance.NAWABI2_dLGN_No1",	"Abundance.NAWABI2_dLGN_No_2",
                    "Abundance.NAWABI2_dLGN_No_3", "Abundance.NAWABI2_dLGN_No_4",
                    "Abundance.NAWABI2_dLGN_Crush_1",	"Abundance.NAWABI2_dLGN_Crush_2",
                    "Abundance.NAWABI2_dLGN_Crush_3",	"Abundance.NAWABI2_dLGN_Crush_4")
  } else if (target == "chiasma") {
    samplesname = c("Abundance.BELIN1_chiasmaNo1", "Abundance.BELIN1_chiasmaNo2",
                    "Abundance.BELIN1_chiasmaNo3", "Abundance.BELIN1_chiasmaNo4",
                    "Abundance.BELIN1_chiasmaCrush1", "Abundance.BELIN1_chiasmaCrush2",
                    "Abundance.BELIN1_chiasmaCrush3", "Abundance.BELIN1_chiasmaCrush4")
  } else if (target == "SC") {
    samplesname = c("Abundance.NAWABI2_SC_No1",	"Abundance.NAWABI2_SC_No_2",
                    "Abundance.NAWABI2_SC_No_3","Abundance.NAWABI2_SC_No_4",
                    "Abundance.NAWABI2_SC_Crush_1", "Abundance.NAWABI2_SC_Crush_2",
                    "Abundance.NAWABI2_SC_Crush_3", "Abundance.NAWABI2_SC_Crush_4")
  } else if (target == "SCN") {
    samplesname = c("Abundance.BELIN1_SCNNo1", "Abundance.BELIN1_SCNNo2",
                    "Abundance.BELIN1_SCNNo3", "Abundance.BELIN1_SCNNo4",
                    "Abundance.BELIN1_SCNCrush1","Abundance.BELIN1_SCNCrush2",
                    "Abundance.BELIN1_SCNCrush3", "Abundance.BELIN1_SCNCrush4")
  } else if (target == "vLGN") {
    samplesname = c("Abundance.NAWAIBI2_vLGN_No1", "Abundance.NAWAIBI2_vLGN_No2",
                    "Abundance.NAWABI2_vLGN_No_3", "Abundance.NAWABI2_vLGN_No_4",
                    "Abundance.NAWAIBI2_vLGN_crush1", "Abundance.NAWAIBI2_vLGN_crush2",
                    "Abundance.NAWABI2_vLGN_Crush_3", "Abundance.NAWABI2_vLGN_Crush_4")
  }
  
  mask_common_prots = grepl(", ", allproteins)
  common_prots = allproteins[mask_common_prots]
  all_common_prots = unlist(sapply(common_prots, strsplit, ", ")) # Get prots from shared peptides
  alluniq_prots = unique(c(all_common_prots, allproteins[!mask_common_prots]))
  
  pep.data = matrix(t(peptide_data[, samplesname]), nrow = length(samplesname))
  colnames(pep.data) = alluniq_pep
  rownames(pep.data) = samplesname
  
  adj = matrix(F, ncol = length(alluniq_prots), nrow = length(alluniq_pep))
  colnames(adj) = alluniq_prots
  rownames(adj) = alluniq_pep
  
  for (i in 1:length(alluniq_prots)) {
    adj[, i] = grepl(paste(alluniq_prots[i], "($|, )", sep = ""), allproteins)
  }
  
  # Remove contaminants
  mask_con_prot = grepl("#C#", alluniq_prots)
  adj = adj[, !mask_con_prot]
  mask_con_pep = rowSums(adj) == 0
  pep.data = pep.data[, !mask_con_pep]
  adj = adj[!mask_con_pep, ]
  
  pep.data[pep.data == 0] = NA
  
  processed.data = list(peptides_ab = pep.data, adj = adj)
  
  if (!is.null(n.na.max)) {
    processed.data = remove_NA_pep_capizzi(processed.data, n.na.max)
  }
  pep.data = processed.data$peptides_ab

  if (is.null(norm.method)) {
    pep.data = log2(pep.data)
  }
  else if (norm.method == "vsn") {
    pep.data = t(vsn::justvsn(t(pep.data)))
  }
  else if (norm.method == "median") {
    pep.data = log2(pep.data)
    sample_medians = apply(pep.data, 1, median,  na.rm = T)
    mean_sample_median = mean(sample_medians)
    pep.data = pep.data + mean_sample_median - sample_medians
  }
  processed.data$peptides_ab = pep.data
  
  if (rm.nested.pg) {
    idx.emb.prots = get_indexes_embedded_prots(processed.data$adj)
    processed.data = rm_pg_from_idx_merge_pg(processed.data, idx.emb.prots)
  }
  
  return(processed.data)
}
