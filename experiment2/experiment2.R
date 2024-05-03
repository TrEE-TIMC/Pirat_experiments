library(foreach) # CRAN
library(doParallel) # CRAN
source("../data_load/capizzi_data.R")
source("../utils.R")

# Imputation methods
library(Pirat) # GitHub
library(msImpute) # Bioconductor
library(MsCoreUtils) # Bioconductor
library(missForest) # CRAN
library(impute) # Bioconductor
library(imputeLCMD) # CRAN
library(pcaMethods) # Bioconductor
library(tidyverse) # CRAN
library(SeqKnn) # CRAN
library(GMSimpute) # CRAN
library(rrcovNA) # CRAN
source('../trKNN/Imput_funcs.R')

# Create repo of all results
dir.create("res")

# Choose the dataset of the experiment for all datasets
#
# Load Capizzi data
# data.name = "Capizzi2022"
# pseudo.mv.rate = 0.16
# npcs = 2
# pep.data.comp = readRDS("../data/capizzi_comp.rds")
# folderexp = file.path("res", data.name) 

# Load Vilallongue SCN data
data.name = "Vilallongue2022_SCN"
pseudo.mv.rate = 0.14
npcs = 2
pep.data.comp = readRDS("../data/vilallongue_scn_ion_vsn.rds")
folderexp = file.path("res", data.name)

# Load Vilallongue SC data
# data.name = "Vilallongue2022_SC"
# pseudo.mv.rate = 0.14
# npcs = 2
# pep.data.comp = readRDS("../data/vilallongue_sc_ion_vsn.rds")
# folderexp = file.path("res", data.name)

n.cores = 2 # Number of cores for imputation in parallel.
seednums = 0:4 + 543210
mnar.mv.rates = seq(0, 1, 0.25) # pseudo MVs rates

impmethods = c("DATA", "Pirat", "Pirat_degenerated", "MinProb", "QRILC", "SeqKNN", "ImpSeq", "BPCA")


nsamples = nrow(pep.data.comp$peptides_ab)
npeps = ncol(pep.data.comp$peptides_ab)
tol.pseudo.na = nsamples - 1 # Peptides with tol.pseudo.na or more MVs are filtered out
mncar.fold = paste("MNAR", mnar.mv.rates)

# Get interpolated curve of Probit mean vs given MNAR proportion
sd.t = sd(pep.data.comp$peptides_ab, na.rm = T)/2 # sd of the probit mechanism
expec_sum_m = Vectorize(
  function(mu, xx, sd) {return(1 - mean(pnorm(xx, mu, sd), na.rm = T))},
  vectorize.args = "mu")
qlow = quantile(pep.data.comp$peptides_ab, probs = pseudo.mv.rate*0.001, na.rm = T)
qhigh = quantile(pep.data.comp$peptides_ab, probs = pseudo.mv.rate, na.rm = T)
q_vec = seq(qlow, qhigh, length.out = 100)
expec_sum_m_vec = expec_sum_m(q_vec, pep.data.comp$peptides_ab, sd.t)
approx_q = approxfun(expec_sum_m_vec, q_vec)

# Make clusters if necessary
if (n.cores > 1) {
  cur.cluster = makeCluster(n.cores, type = "FORK", outfile = "")
  print(cur.cluster)
  registerDoParallel(cl = cur.cluster)
  getDoParRegistered()
  getDoParWorkers()
} else {
  registerDoSEQ()
}

n.seeds = length(seednums)
n.mnars = length(mnar.mv.rates)

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
        q = approx_q(mnar.mv.rate*pseudo.mv.rate)
      }
      tt = matrix(rnorm(nsamples*npeps, q, sd.t), nrow = nsamples,
                  ncol = npeps)
      pseudo.data = pep.data.comp
      under.t = matrix(F, nrow = nsamples,
                       ncol = npeps)
      under.t[pep.data.comp$peptides_ab < tt] = T
      mnar.mask = under.t
      pseudo.mask = mnar.mask
      boolvals = as.logical(rbinom(sum(!mnar.mask), 1,
                                   pseudo.mv.rate*(1 - mnar.mv.rate)
                                   / (1 - pseudo.mv.rate*mnar.mv.rate)))
      pseudo.mask[!mnar.mask] = boolvals
      pseudo.data$peptides_ab[pseudo.mask] = NA
      # print(mean(pseudo.mask))
  
      # plot2hists(pep.data.comp$peptides_ab, pseudo.data$peptides_ab, "Observed values",
                 # "Missing Values", "MCAR scenario", "Peptides abundances", freq = T)
  
      # Remove peptide w too many NAs from pseudo dataset in both control and pseudo datasets
      idxpep.fullna = which(colSums(is.na(pseudo.data$peptides_ab)) >= tol.pseudo.na)
      comp.rmfullna = remove_pep_from_idx(pep.data.comp, idxpep.fullna)
      pseudo.data = remove_NA_pep_capizzi(pseudo.data, tol.pseudo.na)
      pseudo.data[["comp_pep_abs"]] = comp.rmfullna$peptides_ab
  
      # print(paste("NA proportion", mean(is.na(pseudo.data$peptides_ab))))
      # print(paste("Dim adj", paste(dim(pseudo.data$adj))))
      # print(paste("N shared pep", sum(rowSums(pseudo.data$adj) >= 2)))
      # print(paste("N pseudo NA", sum(is.na(pseudo.data$peptides_ab) & !is.na(pseudo.data$comp_pep_abs))))
      
      if (method == "DATA") {
        idx.emb.prots = get_indexes_embedded_prots(pseudo.data$adj)
        pseudo.data = rm_pg_from_idx_merge_pg(pseudo.data, idx.emb.prots)
        saveRDS(pseudo.data, file = file.path(path2saveRDS, "DATA.rds"))
      }
      
      # Pirat
      if (method == "Pirat") {
        start_time <- Sys.time()
        res.mle_mnar = my_pipeline_llkimpute(pseudo.data, 
                                             pep.ab.comp=pseudo.data$comp_pep_abs)
        res.mle_mnar = res.mle_mnar$data.imputed
        end_time = Sys.time()
        saveRDS(res.mle_mnar, file = file.path(path2saveRDS, "Pirat.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "Pirat_time.rds"))
      }
      
      # Pirat degenerated
      if (method == "Pirat_degenerated") {
        start_time <- Sys.time()
        res.mle_mnar = my_pipeline_llkimpute(pseudo.data, 
                                             pep.ab.comp=pseudo.data$comp_pep_abs,
                                             degenerated = T)
        res.mle_mnar = res.mle_mnar$data.imputed
        end_time = Sys.time()
        saveRDS(res.mle_mnar, file = file.path(path2saveRDS, "Pirat_degenerated.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "Pirat_degenerated_time.rds"))
      }
      
      
      # KNN features
      if (method == "KNN") {
        library("impute")
        start_time <- Sys.time()
        res.knn = impute.knn(t(pseudo.data$peptides_ab), k = 10, rowmax = 1, colmax = 1)
        res.knn$data = t(res.knn$data)
        saveRDS(res.knn$data, file = file.path(path2saveRDS, "KNN.rds"))
        end_time = Sys.time()
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "KNN_time.rds"))
      }
      
      # BPCA
      if (method == "BPCA") {
        start_time <- Sys.time()
        res.bpca = pca(t(pseudo.data$peptides_ab), method = "bpca", nPcs = npcs)
        end_time = Sys.time()
        res.bpca = t(res.bpca@completeObs)
        saveRDS(res.bpca, file = file.path(path2saveRDS, "BPCA.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        print(difference)
        saveRDS(difference, file.path(path2saveRDS, "BPCA_time.rds"))
      }
      
      # MissForest
      if (method == "MissForest") {
        start_time <- Sys.time()
        res.misfor = missForest(pseudo.data$peptides_ab)
        saveRDS(res.misfor$ximp, file = file.path(path2saveRDS, "MissForest.rds"))
        end_time = Sys.time()
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "MissForest_time.rds"))
      }
      
      # MLE
      if (method == "MLE") {
        start_time <- Sys.time()
        res.mleclass = t(impute_mle(t(pseudo.data$peptides_ab)))
        end_time = Sys.time()
        saveRDS(res.mleclass, file = file.path(path2saveRDS, "MLE.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "MLE_time.rds"))
      }
      
      # QRILC # Need at least 2 values
      if (method == "QRILC") {
        start_time <- Sys.time()
        res.qrilc = impute.QRILC(t(pseudo.data$peptides_ab))
        res.qrilc = t(res.qrilc[[1]])
        end_time = Sys.time()
        saveRDS(res.qrilc, file = file.path(path2saveRDS, "QRILC.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        print(difference)
        saveRDS(difference, file.path(path2saveRDS, "QRILC_time.rds"))
      }
      
      # MinProb
      if (method == "MinProb") {
        start_time <- Sys.time()
        res.minprob = t(impute.MinProb(t(pseudo.data$peptides_ab)))
        end_time = Sys.time()
        saveRDS(res.minprob, file = file.path(path2saveRDS, "MinProb.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        print(difference)
        saveRDS(difference, file.path(path2saveRDS, "MinProb_time.rds"))
      }
      
      # SVD
      if (method == "SVD") { 
        start_time <- Sys.time()
        res.svd = pca(t(pseudo.data$peptides_ab), nPcs = npcs, method = "svdImpute")
        end_time = Sys.time()
        res.svd = t(res.svd@completeObs)
        saveRDS(res.svd, file = file.path(path2saveRDS, "SVD.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "SVD_time.rds"))
      }
      
      # LLS
      if (method == "LLS") { 
        start_time <- Sys.time()
        res.lls = llsImpute(pseudo.data$peptides_ab)
        end_time = Sys.time()
        res.lls = res.lls@completeObs
        saveRDS(res.lls, file = file.path(path2saveRDS, "LLS.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "LLS_time.rds"))
      }
      
      # trKNN
      if (method == "trKNN") { 
        start_time <- Sys.time()
        res.trknn <- t(imputeKNN(
          t(pseudo.data$peptides_ab[, !colSums(is.na(pseudo.data$peptides_ab)) >= 3]), k = 10, distance = 'truncation'
        ))
        res.trknn = put_imputed_in_or(res.trknn, pseudo.data$peptides_ab)
        end_time = Sys.time()
        saveRDS(res.trknn, file = file.path(path2saveRDS, "trKNN.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "trKNN_time.rds"))
      }

      # Seq-KNN
      if (method == "SeqKNN") { 
        start_time <- Sys.time()
        res.seqknn <- t(SeqKNN(t(pseudo.data$peptides_ab), k = 10))
        end_time = Sys.time()
        saveRDS(res.seqknn, file = file.path(path2saveRDS, "SeqKNN.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "SeqKNN_time.rds"))
      }

      # GMS
      if (method == "GMS") { 
        start_time <- Sys.time()
        res.gms <- t(GMS.Lasso(
          t(
            pseudo.data$peptides_ab[, colSums(!is.na(pseudo.data$peptides_ab)) >= 3]
          ),nfolds=3,log.scale=FALSE,TS.Lasso=TRUE))
        res.gms = put_imputed_in_or(res.gms, pseudo.data$peptides_ab)
        end_time = Sys.time()
        saveRDS(res.gms, file = file.path(path2saveRDS, "GMS.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "GMS_time.rds"))
      }
      
      # Impseq
      if (method == "ImpSeq") { 
        start_time <- Sys.time()
        res.impseq <- t(impSeq(t(pseudo.data$peptides_ab)))
        end_time = Sys.time()
        saveRDS(res.impseq, file = file.path(path2saveRDS, "ImpSeq.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "ImpSeq_time.rds"))
      }

      # ImpSeqRob
      if (method == "ImpSeqRob") { 
        start_time <- Sys.time()
        res.impseqrob <- t(impSeqRob(t(pseudo.data$peptides_ab))$x)
        end_time = Sys.time()
        saveRDS(res.impseqrob, file = file.path(path2saveRDS, "ImpSeqRob.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "ImpSeqRob_time.rds"))
      }
      print(paste(method, "DONE !"))
      NULL
    }
overall_end <- Sys.time()
if (n.cores > 1) {
  stopCluster(cl = cur.cluster)
}
print(overall_end - overall_start)

