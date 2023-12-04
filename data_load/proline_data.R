source("../data_loader/capizzi_data.R")

get.proline.data = function(path, n.na.max = NULL, rm.nested.pg = F) {
  
  peptide_data = read.csv(path, sep = ";", dec = ",")

  alluniq_prots = unique(peptide_data$accession)
  pep_modif = paste(peptide_data$sequence, peptide_data$modifications)
  alluniq_pep = unique(pep_modif)
  
  samplesname = c("raw_abundance_10amol_1", "raw_abundance_10amol_2",
                  "raw_abundance_10amol_3", "raw_abundance_10amol_4",
                  "raw_abundance_50amol_1", "raw_abundance_50amol_2",
                  "raw_abundance_50amol_3", "raw_abundance_50amol_4",
                  "raw_abundance_100amol_1", "raw_abundance_100amol_2",
                  "raw_abundance_100amol_3", "raw_abundance_100amol_4",
                  "raw_abundance_250amol_1", "raw_abundance_250amol_2",
                  "raw_abundance_250amol_3", "raw_abundance_250amol_4",
                  "raw_abundance_500amol_1", "raw_abundance_500amol_2",
                  "raw_abundance_500amol_3", "raw_abundance_500amol_4",
                  "raw_abundance_1fmol_1", "raw_abundance_1fmol_2",
                  "raw_abundance_1fmol_3", "raw_abundance_1fmol_4",
                  "raw_abundance_5fmol_1", "raw_abundance_5fmol_2",
                  "raw_abundance_5fmol_3", "raw_abundance_5fmol_4",
                  "raw_abundance_10fmol_1", "raw_abundance_10fmol_2",
                  "raw_abundance_10fmol_3", "raw_abundance_10fmol_4",
                  "raw_abundance_25fmol_1", "raw_abundance_25fmol_2",
                  "raw_abundance_25fmol_3", "raw_abundance_25fmol_4",
                  "raw_abundance_50fmol_1", "raw_abundance_50fmol_2",
                  "raw_abundance_50fmol_3", "raw_abundance_50fmol_4"
  )
  
  pep.data = matrix(NA, ncol = length(alluniq_pep), nrow = length(samplesname))
  colnames(pep.data) = alluniq_pep
  rownames(pep.data) = samplesname
  
  adj = matrix(0, ncol = length(alluniq_prots), nrow = length(alluniq_pep))
  colnames(adj) = alluniq_prots
  rownames(adj) = alluniq_pep
  
  for (i in 1:nrow(peptide_data)) {
    idx_uniq_pep = which(pep_modif[i] == alluniq_pep)
    idx_uniq_prot = which(peptide_data$accession[i] == alluniq_prots)
    pep.data[, idx_uniq_pep] = t(peptide_data[i, samplesname])
    adj[idx_uniq_pep, idx_uniq_prot] = 1
  }
  
  # Remove contaminants
  mask_con_prot = grepl("#C#", alluniq_prots)
  adj = adj[, !mask_con_prot]
  mask_con_pep = rowSums(adj) == 0
  pep.data = pep.data[, !mask_con_pep]
  adj = adj[!mask_con_pep, ]
  
  mask_human_prot = grepl("_HUMAN_UPS", colnames(adj))
  mask_human_pep = rowSums(adj[, mask_human_prot]) != 0
  
  processed.data = list(peptides_ab = log2(pep.data), adj = adj, 
                        mask_prot_diff = mask_human_prot, mask_pep_diff = mask_human_pep)

  if (!is.null(n.na.max)) {
    processed.data = remove_NA_pep_maxquant(processed.data, n.na.max)
  }
  
  if (rm.nested.pg) {
    idx.emb.prots = get_indexes_embedded_prots(processed.data$adj)
    processed.data = rm_pg_from_idx_merge_pg(processed.data, idx.emb.prots)
  }
  
  return(processed.data)
}
