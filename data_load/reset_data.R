# Get proteogenomic data from RESET or MouseColon. Build PG according to Protein belonging
get.reset.data = function(path, sel_pep, sel_rna) {
  
  new.table <- read.csv(path, sep="\t")
  
  new.table = cbind(new.table[,1:6], new.table[, sel_pep], new.table[, sel_rna]) # select columns of interest
  new.table[,sel_pep][is.na(new.table[,sel_pep])] = 0 
  
  # Les NA proviennent des vecteurs ins?r?s dans la bact?rie, ou des contaminants.
  # Les 0 proviennent d'une vraie abscence de transcrit
  allgenes = new.table$gene_id
  mask.gene.sel = !is.na(allgenes)
  rows.gene.sel = which(mask.gene.sel)
  allgenes.uniq = unique(allgenes[mask.gene.sel])
  allprot.uniq = unique(new.table$prot_id[mask.gene.sel])
  allpeps.uniq = unique(new.table$Sequence_PTM_Charge[mask.gene.sel])
  
  pep.data = matrix(0, nrow = length(sel_pep), ncol = length(allpeps.uniq))
  colnames(pep.data) = allpeps.uniq
  rownames(pep.data) = sel_pep
  
  rna.data = matrix(0, nrow = length(sel_rna), ncol = length(allgenes.uniq))
  colnames(rna.data) = allgenes.uniq
  rownames(rna.data) = sel_rna
  
  adj_pep = matrix(0, nrow = length(allpeps.uniq), ncol = length(allprot.uniq))
  colnames(adj_pep) = allprot.uniq
  rownames(adj_pep) = allpeps.uniq
  
  adj_rna = matrix(0, nrow = length(allgenes.uniq), ncol = length(allprot.uniq))
  colnames(adj_rna) = allprot.uniq
  rownames(adj_rna) = allgenes.uniq
  
  charges = rep(0, length(rows.gene.sel))
  modifs = rep("", length(rows.gene.sel))
  
  for (i in rows.gene.sel) {
    ion.cur = new.table$Sequence_PTM_Charge[i]
    ion.ab.cur = as.numeric(new.table[i, sel_pep])
    
    gene.cur = new.table$gene_id[i]
    prot.cur = new.table$prot_id[i]
    rna.ab.cur = as.numeric(new.table[i, sel_rna])
    
    charge.cur = new.table$master_quant_peptide_ion_charge[i]
    modif.cur = new.table$modifications[i]
    
    mask.cur.ion = which(allpeps.uniq == ion.cur)
    pep.ab2change = pep.data[,mask.cur.ion]
    # stopifnot(all(pep.ab2change == 0) | all(pep.ab2change == ion.ab.cur))
    if (!(all(pep.ab2change == 0) | all(pep.ab2change == ion.ab.cur))) {
      if (max(abs(pep.ab2change - ion.ab.cur)) >= 1e-3) {
        print(pep.ab2change)
        print(ion.ab.cur)
        break
      }
    }
    pep.data[, mask.cur.ion] = ion.ab.cur
    
    
    mask.cur.gene = which(allgenes.uniq == gene.cur)
    mask.cur.prot = which(allprot.uniq == prot.cur)
    rna.ab2change = rna.data[,mask.cur.gene]
    if (!(all(rna.ab2change == 0) | all(rna.ab2change == rna.ab.cur))) {
      if (max(abs(rna.ab2change - rna.ab.cur)) >= 1e-3) {
        print(rna.ab2change)
        print(rna.ab.cur)
        break
      }
    }
    rna.data[,mask.cur.gene] = rna.ab.cur
    
    adj_pep[mask.cur.ion, mask.cur.prot] = 1
    adj_rna[mask.cur.gene, mask.cur.prot] = 1
    
    charges[mask.cur.ion] = charge.cur
    
    modifs[mask.cur.ion] = modif.cur
  }
  
  pep.data[pep.data == 0] = NA
  # rna.data[rna.data == 0] = runif(sum(rna.data == 0), 0.5, 1.5)
  
  return(list(peptides_ab = pep.data, rnas_ab = log2(rna.data+1), adj = adj_pep,
              adj_rna_pg = adj_rna, charges = charges, modifs = modifs))
}


# Remove peptides with number of MVs higher or equal to t
remove_NA_pep_reset <- function(l_pep_rna, t, percentage=F) {
  nsample = nrow(l_pep_rna$rnas_ab)
  if (percentage) {
    t = ceiling(nsample*t)
  }
  mv_count_pep = colSums(is.na(l_pep_rna$peptides_ab))
  i_pep_rm = which(mv_count_pep >= t)
  l_pep_rna = remove_pep_from_idx(l_pep_rna, i_pep_rm)
  return(l_pep_rna)
}

#' Removes 
#'
#' @param l_pep_rna 
#' @param cond_idx 
#' @param n_cond_min 
#' @param offset_log 
#'
#' @return
#' @export
#'
#' @examples
remove_Os_rna_reset <- function(l_pep_rna, cond_idx, n_cond_min=2,
                                offset_log = 1) {
  rnas_ab = 2^l_pep_rna$rnas_ab - offset_log
  count_present = matrix(F, length(unique(cond_idx)), ncol(rnas_ab))
  for (cond in unique(cond_idx)) {
    count_present[cond, ] = colSums(rnas_ab[cond_idx == cond, , drop = F] > 0) >= 1
  }
  idx_rna_2_rm = which(colSums(count_present) < n_cond_min)
  if (length(idx_rna_2_rm) != 0) {
    l_pep_rna$adj_rna_pg = l_pep_rna$adj_rna_pg[-idx_rna_2_rm, ]
    l_pep_rna$rnas_ab = l_pep_rna$rnas_ab[, -idx_rna_2_rm]
    return(l_pep_rna)
  }
  else {
    return(l_pep_rna)
  }

}



# Remove shared peptides
remove_shared_pep <- function(l_pep_rna) {
  nsample = nrow(l_pep_rna$rnas_ab)
  i_pep_rm = which(rowSums(l_pep_rna$adj) >= 2)
  if (length(i_pep_rm) != 0) {
    new_adj = l_pep_rna$adj[-i_pep_rm,]
    i_rna_rm = which(colSums(new_adj) == 0)
    new_pep = l_pep_rna$peptides_ab[,-i_pep_rm]
    new_z = l_pep_rna$charges[-i_pep_rm]
    new_modifs = l_pep_rna$modifs[-i_pep_rm]
    if (length(i_rna_rm) != 0) {
      new_rna = l_pep_rna$rnas_ab[,-i_rna_rm]
      new_adj = new_adj[,-i_rna_rm]
      return(list("peptides_ab" = new_pep, "rnas_ab" = new_rna, "adj" = new_adj,
                  "charges" = new_z, "modifs" = new_modifs))
    }
    else {
      return(list("peptides_ab" = new_pep, "rnas_ab" = l_pep_rna$rnas_ab,
                  "adj" = new_adj, "charges" = new_z, "modifs" = new_modifs))
    }
    
  }
  else {
    return(l_pep_rna)
  }
}