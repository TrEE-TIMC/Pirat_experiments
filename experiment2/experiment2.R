setwd("~/PhD_project/2020-proteomics-transcriptomics/experiment_script/")
source("../data_loader/capizzi_data.R")
source("../data_loader/vilallongue_data.R")
library(reticulate)
library(foreach)
library(doParallel)

put_imputed_in_or = function(imputed.ds, or.ds) {
  imputed.names = colnames(imputed.ds)
  or.names = colnames(or.ds)
  is.imputed = or.names %in% imputed.names
  or.ds[, is.imputed] = imputed.ds
  return(or.ds)
}


# Load Capizzi data
data.name = "Capizzi2022"
pseudo.mv.rate = 0.16
npcs = 2
pep.data.comp = readRDS("../processed_data/capizzi_comp.rds")
folderexp = "pseudo_na_percent_"


# Load Vilallongue data
data.name = "Vilallongue2022"
pseudo.mv.rate = 0.14
npcs = 2
pep.data.comp = readRDS("../processed_data/vilallongue_scn_ion_vsn.rds")
folderexp = "pseudo_SCN_ion_vsn_na_percent_"

n.cores = 2 # Number of cores for imputation in parallel.
seednums = 0:4 + 543210
mnar.mv.rates = c(0, 0.25, 0.5, 0.75, 1)

impmethods = c("KNN", "GMS", "ImpSeqRob", "Pirat", "MissForest", "MinProb",
               "QRILC", "SVD", "LLS", "trKNN", "SeqKNN", "ImpSeq", "BPCA", 
               "MLE", "msImpute_mar", "msImpute_mnar")


nsamples = nrow(pep.data.comp$peptides_ab)
npeps = ncol(pep.data.comp$peptides_ab)

# Add pseudo missing values
tol.pseudo.na = nsamples - 1
mncar.fold = paste("MNAR", mnar.mv.rates)
sd.t = sd(pep.data.comp$peptides_ab, na.rm = T)/2
# Get approximated function of censoring mechanism mean vs given MNAR proportion
expec_sum_m = Vectorize(
  function(mu, xx, sd) {return(1 - mean(pnorm(xx, mu, sd), na.rm = T))},
  vectorize.args = "mu")
qlow = quantile(pep.data.comp$peptides_ab, probs = pseudo.mv.rate*0.001, na.rm = T)
qhigh = quantile(pep.data.comp$peptides_ab, probs = pseudo.mv.rate, na.rm = T)
q_vec = seq(qlow, qhigh, length.out = 100)
expec_sum_m_vec = expec_sum_m(q_vec, pep.data.comp$peptides_ab, sd.t)
approx_q = approxfun(expec_sum_m_vec, q_vec)


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

path2saveexp = file.path("..", "experiments", data.name,
                         paste(folderexp, pseudo.mv.rate*100,
                               sep=""))
dir.create(path2saveexp)

overall_start <- Sys.time()
foreach(seednum = seednums) %:%
  foreach(idx.mnar.mv.rate = 1:n.mnars) %:%
    foreach(method = impmethods) %dopar% {
      print("###################")
      print(paste(seednum, idx.mnar.mv.rate, method))
      print("###################")
      mnar.mv.rate = mnar.mv.rates[idx.mnar.mv.rate]
      path2saveRDS = file.path(path2saveexp, mncar.fold[idx.mnar.mv.rate])
      dir.create(path2saveRDS)
      path2saveRDS = file.path(path2saveRDS, seednum)
      dir.create(path2saveRDS)
      set.seed(seednum)
      if (mnar.mv.rate == 0) {
        q = -1000
      } else {
        q = approx_q(mnar.mv.rate*pseudo.mv.rate) #as.logical(rbinom(sum(under.t), 1, mnar.mv.rate))
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
      
      # Our method
      if (method == "MLEMNAR") {
        source("../method/pipeline_impute.R")
        start_time <- Sys.time()
        res.mle_mnar = pipeline_llkimpute(pseudo.data, 
                                          pep.ab.comp=pseudo.data$comp_pep_abs, )
        end_time = Sys.time()
        saveRDS(res.mle_mnar, file = file.path(path2saveRDS, "MLEMNAR_2.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "MLEMNAR_2_time.rds"))
      }
      
      # Our method transpose
      if (method == "MLEMNAR_transpose") {
        source("../method/pipeline_impute.R")
        start_time <- Sys.time()
        res.mle_mnar = pipeline_llkimpute(pseudo.data, 
                                          pep.ab.comp=pseudo.data$comp_pep_abs,
                                          nu_factor=2, protidxs = NULL,
                                          transpose = T)
        end_time = Sys.time()
        saveRDS(res.mle_mnar, file = file.path(path2saveRDS, "MLEMNAR_transpose.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "MLEMNAR_transpose_time.rds"))
      }
      
      # Our method transpose
      if (method == "MLEMNAR_degenerated") {
        source("../method/pipeline_impute.R")
        start_time <- Sys.time()
        res.mle_mnar = pipeline_llkimpute(pseudo.data, 
                                          pep.ab.comp=pseudo.data$comp_pep_abs,
                                          degenerated = T)
        end_time = Sys.time()
        saveRDS(res.mle_mnar, file = file.path(path2saveRDS, "MLEMNAR_degenerated.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "MLEMNAR_degenerated_time.rds"))
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
        library("pcaMethods")
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
        library("missForest")
        start_time <- Sys.time()
        res.misfor = missForest(pseudo.data$peptides_ab)
        saveRDS(res.misfor$ximp, file = file.path(path2saveRDS, "MissForest.rds"))
        end_time = Sys.time()
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "MissForest_time.rds"))
      }
      
      # MLE classic
      if (method == "MLECLASS") {
        library("MsCoreUtils")
        start_time <- Sys.time()
        res.mleclass = t(impute_mle(t(pseudo.data$peptides_ab)))
        end_time = Sys.time()
        saveRDS(res.mleclass, file = file.path(path2saveRDS, "MLECLASS.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "MLECLASS_time.rds"))
      }
      
      # QRILC # Need at least 2 values
      if (method == "QRILC") {
        library("imputeLCMD")
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
        library("imputeLCMD")
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
        library(pcaMethods)
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
        library(pcaMethods)
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
        source('../Imput_funcs.R')
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
      if (method == "seqKNN") { 
        library(SeqKnn)
        start_time <- Sys.time()
        res.seqknn <- t(SeqKNN(t(pseudo.data$peptides_ab), k = 10))
        end_time = Sys.time()
        saveRDS(res.seqknn, file = file.path(path2saveRDS, "seqKNN.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "seqKNN_time.rds"))
      }

      # GMS
      if (method == "GMS") { 
        library(GMSimpute)
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
        library(rrcovNA)
        start_time <- Sys.time()
        res.impseq <- t(impSeq(t(pseudo.data$peptides_ab)))
        end_time = Sys.time()
        saveRDS(res.impseq, file = file.path(path2saveRDS, "ImpSeq.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "ImpSeq_time.rds"))
      }

      # ImpSeqRob
      if (method == "ImpSeqRob") { 
        library(rrcovNA)
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

