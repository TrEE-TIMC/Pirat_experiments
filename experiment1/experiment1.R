library(Pirat)
library(msImpute)
library(MsCoreUtils)
library(missForest)
library(impute)
library(imputeLCMD)
library(pcaMethods)
library(tidyverse)
library(SeqKnn)
library(GMSimpute)
library(rrcovNA)

source('../trKNN/Imput_funcs.r')


put_imputed_in_or = function(imputed.ds, or.ds) {
  imputed.names = colnames(imputed.ds)
  or.names = colnames(or.ds)
  is.imputed = or.names %in% imputed.names
  or.ds[, is.imputed] = imputed.ds
  return(or.ds)
}


# Choose the dataset to impute by uncommenting following 3 blocks

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
start_time <- Sys.time()
res.mle_mnar = pipeline_llkimpute(mq.data.comp, pep.ab.comp=NULL, nu_factor=2)
res.mle_mnar = res.mle_mnar$data.imputed
end_time = Sys.time()
saveRDS(res.mle_mnar, file = file.path(res.fold.name, "Pirat.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "Pirat_time.rds"))

# MsImpute MNAR
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
start_time <- Sys.time()
res.mleclass = impute_mle(mq.data.comp$peptides_ab)
end_time = Sys.time()
saveRDS(res.mleclass, file = file.path(res.fold.name, "MLE.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "MLE_time.rds"))


# MsImpute
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
start_time <- Sys.time()
res.misfor = missForest(mq.data.comp$peptides_ab)
end_time = Sys.time()
saveRDS(res.misfor$ximp, file = file.path(res.fold.name, "MissForest.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "MissForest_time.rds"))

# KNN
start_time <- Sys.time()
res.knn = impute.knn(t(mq.data.comp$peptides_ab), k = 10, rowmax = 1, colmax = 1)
end_time = Sys.time()
res.knn = t(res.knn$data)
saveRDS(res.knn, file = file.path(res.fold.name, "KNN.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "KNN_time.rds"))


# QRILC
start_time <- Sys.time()
res.qrilc = impute.QRILC(t(mq.data.comp$peptides_ab))
end_time = Sys.time()
res.qrilc = t(res.qrilc[[1]])
saveRDS(res.qrilc, file = file.path(res.fold.name, "QRILC.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "QRILC_time.rds"))


# MinProb
start_time <- Sys.time()
res.minprob = t(impute.MinProb(t(mq.data.comp$peptides_ab)))
end_time = Sys.time()
saveRDS(res.minprob, file = file.path(res.fold.name, "MinProb.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "MinProb_time.rds"))

# SVD
start_time <- Sys.time()
res.svd = pca(t(mq.data.comp$peptides_ab), nPcs = 2, method = "svdImpute")
end_time = Sys.time()
res.svd = t(res.svd@completeObs)
saveRDS(res.svd, file = file.path(res.fold.name, "SVD.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "SVD_time.rds"))

# LLS
start_time <- Sys.time()
res.lls = llsImpute(mq.data.comp$peptides_ab)
end_time = Sys.time()
res.lls = res.lls@completeObs
saveRDS(res.lls, file = file.path(res.fold.name, "LLS.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "LLS_time.rds"))

# trKNN
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
start_time <- Sys.time()
res.seqknn <- t(SeqKNN(t(mq.data.comp$peptides_ab), k = 10))
end_time = Sys.time()
saveRDS(res.seqknn, file = file.path(res.fold.name, "SeqKNN.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "SeqKNN_time.rds"))


# GMS
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
start_time <- Sys.time()
res.impseq <- t(impSeq(t(mq.data.comp$peptides_ab)))
end_time = Sys.time()
saveRDS(res.impseq, file = file.path(res.fold.name, "ImpSeq.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "ImpSeq_time.rds"))


# ImpSeqRob
start_time <- Sys.time()
res.impseqrob <- t(impSeqRob(t(mq.data.comp$peptides_ab))$x)
end_time = Sys.time()
saveRDS(res.impseqrob, file = file.path(res.fold.name, "ImpSeqRob.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "ImpSeqRob_time.rds"))


# BPCA
start_time <- Sys.time()
res.bpca = pca(t(mq.data.comp$peptides_ab), method = "bpca", nPcs = npcs)
end_time = Sys.time()
res.bpca = t(res.bpca@completeObs)
saveRDS(res.bpca, file = file.path(res.fold.name, "BPCA.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(res.fold.name, "BPCA_time.rds"))
