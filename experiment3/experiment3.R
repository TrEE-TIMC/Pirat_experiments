library(foreach)
library(doParallel)
library(Pirat)
source("../data_load/reset_data.R")
source("../utils.R")

# Create repo of all results
dir.create("res")

# Load RESET data
name.data = "Ropers2021"
cond_idx_rna = 1:18
cond_idx_pep = 1:18
npcs = 6
all.data.comp = readRDS("../data/ropers_comp_pg.rds")


# Load MouseColon data
# name.data = "Habowski2020"
# npcs = 6
# cond_idx_rna = c(1,1,1,2,2,2,2,3,3,3,3,3,4,4,4,5,5,6,6,6,6,6)
# cond_idx_pep = rep(1:6, each=3)
# all.data.comp = readRDS("../data/habowski_comp_pg.rds")

# Experiment design
impmethods = c("DATA", "Pirat", "Pirat-T", "Pirat-S") # Imputation methods
n.cores = 2 # Number of cores for imputation in parallel.
seednums = 0:9 + 543210
mnar.mv.rates = c(0, 1) # Desired MNAR rate among pseudo MVs
mncar.fold = paste("MNAR", mnar.mv.rates)
nsamples = nrow(all.data.comp$peptides_ab)
npeps = ncol(all.data.comp$peptides_ab)
prop.na = mean(is.na(all.data.comp$peptides_ab))
tol.pseudo.na = nsamples - 1 # Peptides with tol.pseudo.na or more MVs are filtered out

# Get interpolated curve of Probit mean vs given MNAR proportion
n.na = prop.na * npeps * nsamples
n.pseudo.na = 0.01 * n.na # Number of MVs equal to 1% of total number of MVs
pseudo.mv.rate = n.pseudo.na / ((1 - prop.na) * npeps * nsamples)
sd.t = sd(all.data.comp$peptides_ab, na.rm = T)/2
expec_sum_m = Vectorize(
  function(mu, xx, sd) {return(1 - mean(pnorm(xx, mu, sd), na.rm = T))},
  vectorize.args = "mu")
qlow = quantile(all.data.comp$peptides_ab, probs = pseudo.mv.rate*0.001, na.rm = T)
qhigh = quantile(all.data.comp$peptides_ab, probs = pseudo.mv.rate, na.rm = T)
q_vec = seq(qlow, qhigh, length.out = 100)
expec_sum_m_vec = expec_sum_m(q_vec, all.data.comp$peptides_ab, sd.t)
approx_q = approxfun(expec_sum_m_vec, q_vec)

# Make clusters if necessary
if (n.cores > 1) {
  cur.cluster = makeCluster(n.cores, type = "PSOCK", outfile = "")
  print(cur.cluster)
  registerDoParallel(cl = cur.cluster)
  getDoParRegistered()
  getDoParWorkers()
} else {
  registerDoSEQ()
}

n.seeds = length(seednums)
n.mnars = length(mnar.mv.rates)

folderexp = file.path("res", name.data)
dir.create(folderexp)

overall_start <- Sys.time()
foreach(seednum = seednums) %:%
  foreach(idx.mnar.mv.rate = 1:n.mnars) %:%
    foreach(method = impmethods) %dopar% {
      print("###################")
      print(paste(seednum, idx.mnar.mv.rate, method))
      print("###################")
      mnar.mv.rate = mnar.mv.rates[idx.mnar.mv.rate]
      path2saveRDS = file.path(folderexp, mncar.fold[idx.mnar.mv.rate])
      dir.create(path2saveRDS)
      path2saveRDS = file.path(path2saveRDS, seednum)
      dir.create(path2saveRDS)
      set.seed(seednum)
      
      # Add pseudo MVs
      if (mnar.mv.rate == 0) {
        q = -1000
      } else {
        q = approx_q(mnar.mv.rate*pseudo.mv.rate) #as.logical(rbinom(sum(under.t), 1, mnar.mv.rate))
      }
      tt = matrix(rnorm(nsamples*npeps, q, sd.t), nrow = nsamples,
                  ncol = npeps)
      all.pseudo = all.data.comp
      under.t = matrix(F, nrow = nsamples,
                       ncol = npeps)
      under.t[all.data.comp$peptides_ab < tt] = T
      mnar.mask = under.t
      pseudo.mask = mnar.mask
      boolvals = as.logical(rbinom(sum(!mnar.mask), 1,
                                   pseudo.mv.rate*(1 - mnar.mv.rate)
                                   / (1 - pseudo.mv.rate*mnar.mv.rate)))
      pseudo.mask[!mnar.mask] = boolvals
      all.pseudo$peptides_ab[pseudo.mask] = NA
      # print(mean(pseudo.mask))
      
      # plot2hists(all.data.comp$peptides_ab, all.pseudo$peptides_ab, "Wo censoring", "W censoring", "Example", freq = T)
      
      # Remove peptides w too many NAs from pseudo dataset in both control and pseudo datasets
      idxpep.fullna = which(colSums(is.na(all.pseudo$peptides_ab)) >= tol.pseudo.na)
      reset.comp.rmfullna = remove_pep_from_idx(all.data.comp, idxpep.fullna)
      all.pseudo = remove_NA_pep_reset(all.pseudo, tol.pseudo.na)
      all.pseudo[["comp_pep_abs"]] = reset.comp.rmfullna$peptides_ab
      
      # print(paste("NA proportion", mean(is.na(all.pseudo$peptides_ab))))
      # print(paste("Dim adj", paste(dim(all.pseudo$adj))))
      # print(paste("N shared pep", sum(rowSums(all.pseudo$adj) >= 2)))
      # print(paste("N pseudo NA", sum(is.na(all.pseudo$peptides_ab) & !is.na(all.pseudo$comp_pep_abs))))
      
      if (method == "DATA") {
        print("Remove nested prots x2...")
        idx.emb.prots = get_indexes_embedded_prots(all.pseudo$adj)
        all.pseudo = rm_pg_from_idx_merge_pg(all.pseudo, idx.emb.prots)
        saveRDS(all.pseudo, file = file.path(path2saveRDS, "DATA.rds"))
      }
      
      
      # Pirat
      if (method == "Pirat") {
        library(Pirat)
        start_time <- Sys.time()
        res.mle_mnar = pipeline_llkimpute(all.pseudo, 
                                          pep.ab.comp = all.pseudo$comp_pep_abs)
        res.mle_mnar = res.mle_mnar$data.imputed
        end_time = Sys.time()
        saveRDS(res.mle_mnar, file = file.path(path2saveRDS, "Pirat.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "Pirat_time.rds"))
      }
      
      # Pirat-S
      if (method == "Pirat-S") {
        library(Pirat)
        start_time <- Sys.time()
        res.mle_mnar.tpg1 = pipeline_llkimpute(all.pseudo,
                                               pep.ab.comp = all.pseudo$comp_pep_abs,
                                               extension = "S")
        res.mle_mnar.tpg1 = res.mle_mnar.tpg1$data.imputed
        end_time = Sys.time()
        saveRDS(res.mle_mnar.tpg1, file = file.path(path2saveRDS, "Pirat-S.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "Pirat-S.rds"))
      }
      
      # Our method ProteoGenomics
      if (method == "Pirat-T") {
        library(Pirat)
        start_time <- Sys.time()
        res.mle_mnar_pg = pipeline_llkimpute(all.pseudo,
                                             pep.ab.comp = all.pseudo$comp_pep_abs,
                                             extension = "T",
                                             rna.cond.mask = cond_idx_rna,
                                             pep.cond.mask = cond_idx_pep,
                                             max.pg.size.pirat.t = 1)
        res.mle_mnar_pg = res.mle_mnar_pg$data.imputed
        end_time = Sys.time()
        saveRDS(res.mle_mnar_pg, file = file.path(path2saveRDS, "Pirat-T.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "Pirat-T.rds"))
      }
      print(paste(method, "DONE !"))
      NULL
    }
overall_end <- Sys.time()
if (n.cores > 1) {
  stopCluster(cl = cur.cluster)
}
print(overall_end - overall_start)
