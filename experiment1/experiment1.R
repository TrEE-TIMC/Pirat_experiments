# source("../../Pirat/R/pipeline_impute.R")
# source("../../Pirat/R/utils.R")

put_imputed_in_or = function(imputed.ds, or.ds) {
  imputed.names = colnames(imputed.ds)
  or.names = colnames(or.ds)
  is.imputed = or.names %in% imputed.names
  or.ds[, is.imputed] = imputed.ds
  return(or.ds)
}

# Choose the dataset to impute by uncommenting following block

# Load Bouyssié2020 data
#
data.name = "Bouyssie2020"
npcs = 10
mq.data.comp = readRDS(file.path("..", "data", "bouyssie_comp.rds"))
groups = factor(rep(1:10, each = 4))

# Load Cox2014 data
#
# data.name = "Cox2014"
# npcs = 2
# mq.data.comp = readRDS(file.path("..", "processed_data", "cox_comp.rds"))
# groups = factor(rep(1:2, each = 4))

# Load Huang2020 data
#
# data.name = "Huang2020"
# npcs = 5
# mq.data.comp = readRDS(file.path("..", "processed_data", "Huang2020_data_comp.rds"))
# groups = factor(rep(1:5, each = 5))

npeps = ncol(mq.data.comp$peptides_ab)
seednum = 543210
set.seed(seednum)
res.fold.name = "res"
dir.create(res.fold.name)
res.fold.name = file.path(res.fold.name, data.name)
dir.create(res.fold.name)

# Pirat
library(Pirat)
start_time <- Sys.time()
res.mle_mnar = pipeline_llkimpute(mq.data.comp, pep.ab.comp=NULL, nu_factor=2)
end_time = Sys.time()
# res.mle_mnar = readRDS("../../2020-proteomics-transcriptomics/experiments/Proline/543210_w_shared/MLEMNAR.rds")
saveRDS(res.mle_mnar, file = file.path(res.fold.name, "Pirat.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "Pirat_time.rds"))

# MsImpute MNAR
library(msImpute)
start_time <- Sys.time()
res.msimpute_mnar = t(msImpute(
  t(
    mq.data.comp$peptides_ab[, colSums(!is.na(mq.data.comp$peptides_ab)) >= 4]
  ), method = "v2-mnar", group = groups))
res.msimpute_mnar = put_imputed_in_or(res.msimpute_mnar, mq.data.comp$peptides_ab)
end_time = Sys.time()
saveRDS(res.msimpute_mnar, file = file.path(res.fold.name, "msImpute_mnar.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "msImpute_mnar_time.rds"))

# MLE classics
library("MsCoreUtils")
start_time <- Sys.time()
res.mleclass = impute_mle(mq.data.comp$peptides_ab)
end_time = Sys.time()
saveRDS(res.mleclass, file = file.path(res.fold.name, "MLECLASS.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "MLECLASS_time.rds"))


# MsImpute
library(msImpute)
start_time <- Sys.time()
res.msimpute = t(msImpute(
  t(
    mq.data.comp$peptides_ab[, colSums(!is.na(mq.data.comp$peptides_ab)) >= 4]
  ), method = "v2", group = groups))
res.msimpute = put_imputed_in_or(res.msimpute, mq.data.comp$peptides_ab)
end_time = Sys.time()
saveRDS(res.msimpute, file = file.path(res.fold.name, "msImpute.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "msImpute_time.rds"))

# MissForest
library("missForest")
start_time <- Sys.time()
res.misfor = missForest(mq.data.comp$peptides_ab)
end_time = Sys.time()
saveRDS(res.misfor$ximp, file = file.path(res.fold.name, "MissForest.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "MissForest_time.rds"))

# KNN
library("impute")
start_time <- Sys.time()
res.knn = impute.knn(t(mq.data.comp$peptides_ab), k = 10, rowmax = 1, colmax = 1)
end_time = Sys.time()
res.knn = t(res.knn$data)
saveRDS(res.knn, file = file.path(res.fold.name, "KNN.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "KNN_time.rds"))


# QRILC
library("imputeLCMD")
start_time <- Sys.time()
res.qrilc = impute.QRILC(t(mq.data.comp$peptides_ab))
end_time = Sys.time()
res.qrilc = t(res.qrilc[[1]])
saveRDS(res.qrilc, file = file.path(res.fold.name, "QRILC.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "QRILC_time.rds"))


# MinProb
library("imputeLCMD")
start_time <- Sys.time()
res.minprob = t(impute.MinProb(t(mq.data.comp$peptides_ab)))
end_time = Sys.time()
saveRDS(res.minprob, file = file.path(res.fold.name, "MinProb.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "MinProb_time.rds"))

# SVD
library(pcaMethods)
start_time <- Sys.time()
res.svd = pca(t(mq.data.comp$peptides_ab), nPcs = 2, method = "svdImpute")
end_time = Sys.time()
res.svd = t(res.svd@completeObs)
saveRDS(res.svd, file = file.path(res.fold.name, "SVD.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "SVD_time.rds"))

# LLS
library(pcaMethods)
start_time <- Sys.time()
res.lls = llsImpute(mq.data.comp$peptides_ab)
end_time = Sys.time()
res.lls = res.lls@completeObs
saveRDS(res.lls, file = file.path(res.fold.name, "LLS.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "LLS_time.rds"))

# trKNN
source('Imput_funcs.r')
library(tidyverse)
sim_trKNN_wrapper <- function(data) {
  result <- data %>% as.matrix %>% t %>% imputeKNN(., k=10, distance='truncation') %>% t
  return(result)
}
start_time <- Sys.time()
res.trknn <- sim_trKNN_wrapper(
  mq.data.comp$peptides_ab[, !colSums(is.na(mq.data.comp$peptides_ab)) >= 3])
res.trknn = put_imputed_in_or(res.trknn, mq.data.comp$peptides_ab)
end_time = Sys.time()
saveRDS(res.trknn, file = file.path(res.fold.name, "trKNN.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "trKNN_time.rds"))

# Seq-KNN
library(SeqKnn)
start_time <- Sys.time()
res.seqknn <- t(SeqKNN(t(mq.data.comp$peptides_ab), k = 10))
end_time = Sys.time()
saveRDS(res.seqknn, file = file.path(res.fold.name, "seqKNN.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "seqKNN_time.rds"))


# GMS
library(GMSimpute)
start_time <- Sys.time()
res.gms <- t(GMS.Lasso(
  t(
    mq.data.comp$peptides_ab[, colSums(!is.na(mq.data.comp$peptides_ab)) >= 3]
  ),nfolds=3,log.scale=FALSE,TS.Lasso=TRUE))
res.gms = put_imputed_in_or(res.gms, mq.data.comp$peptides_ab)
end_time = Sys.time()
saveRDS(res.gms, file = file.path(res.fold.name, "GMS.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "GMS_time.rds"))


# Impseq
library(rrcovNA)
start_time <- Sys.time()
res.impseq <- t(impSeq(t(mq.data.comp$peptides_ab)))
end_time = Sys.time()
saveRDS(res.impseq, file = file.path(res.fold.name, "ImpSeq.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "ImpSeq_time.rds"))


# ImpSeqRob
library(rrcovNA)
start_time <- Sys.time()
res.impseqrob <- t(impSeqRob(t(mq.data.comp$peptides_ab))$x)
end_time = Sys.time()
saveRDS(res.impseqrob, file = file.path(res.fold.name, "ImpSeqRob.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "ImpSeqRob_time.rds"))


# BPCA
library("pcaMethods")
start_time <- Sys.time()
res.bpca = pca(t(mq.data.comp$peptides_ab), method = "bpca", nPcs = npcs)
end_time = Sys.time()
res.bpca = t(res.bpca@completeObs)
saveRDS(res.bpca, file = file.path(res.fold.name, "BPCA.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "BPCA_time.rds"))
