source("../data_loader/maxquant_data.R")

get.huang.data = function(report_path, n.na.max = NULL, rm.nested.pg = F) {
  dia <- read.delim(report_path, header = TRUE, stringsAsFactors = FALSE)
  
  print("data loaded")
  
  dia$experiment <- paste(dia$R.Condition, dia$R.Replicate,sep="_")
  
  keep2 <- (dia$EG.Qvalue < 0.1)
  
  dia <- dia[keep2,]
  
  dia_ms1 <- tidyr::spread(dia[,c("experiment","EG.PrecursorId", 
                                  "FG.NormalizedMS1PeakArea", 
                                  "PG.ProteinAccessions", "SpikeIn")],
                           key = experiment, value = FG.NormalizedMS1PeakArea)
  allproteins = dia_ms1$PG.ProteinAccessions
  mask_pep_diff = dia_ms1$SpikeIn == "True"
  dia_ms1[dia_ms1==1] <- NA
  
  y_dia_ms1 <- log2(data.matrix(dia_ms1[,-c(1,2,3)]))
  rownames(y_dia_ms1) <- dia_ms1$EG.PrecursorId
  
  y_dda <- y_dia_ms1
  
  y_dda <- t(y_dda)
  
  mask_common_prots = grepl(";", allproteins)
  common_prots = allproteins[mask_common_prots]
  all_common_prots = unlist(sapply(common_prots, strsplit, ";")) # Get prots from shared peptides
  uniq_prots = unique(c(all_common_prots, allproteins[!mask_common_prots]))
  mask_prot_diff = grepl("_UPS", uniq_prots)
  adj = matrix(F, nrow = ncol(y_dda), ncol = length(uniq_prots))
  colnames(adj) = uniq_prots
  rownames(adj) = colnames(y_dda)
  for (i in 1:length(uniq_prots)) {
    adj[, i] = grepl(paste(uniq_prots[i], "($|;)", sep = ""), allproteins)
  }
  
  processed.data = list(peptides_ab = y_dda, adj = adj, 
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




