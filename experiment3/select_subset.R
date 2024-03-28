all.data.comp = readRDS("../data/ropers_comp_pg.rds")
npgs = ncol(all.data.comp$adj)
npgs.sel = 10
idx.pg.sel = sample(npgs, npgs.sel)
idx.sing.pgs = which(colSums(all.data.comp$adj[, idx.pg.sel]) == 1)
idx.pep.sing = which(rowSums(all.data.comp$adj[, idx.pg.sel[idx.sing.pgs]]) >= 1)
idx.pep.sing
all.data.comp$peptides_ab[, idx.pep.sing]
colSums(all.data.comp$adj[, idx.pg.sel])
idx.pep.sel = which(rowSums(all.data.comp$adj[, idx.pg.sel]) >= 1)
peptides_ab = all.data.comp$peptides_ab[, idx.pep.sel]
adj = all.data.comp$adj[idx.pep.sel, idx.pg.sel]
colSums(adj)
idx.rnas.sel = which(rowSums(all.data.comp$adj_rna_pg[, idx.pg.sel]) >= 1)
rnas_ab = all.data.comp$rnas_ab[, idx.rnas.sel]
adj_rna_pg = all.data.comp$adj_rna_pg[idx.rnas.sel, idx.pg.sel]
ropers = list(peptides_ab = peptides_ab, adj = adj, rnas_ab = rnas_ab, adj_rna_pg = adj_rna_pg)
