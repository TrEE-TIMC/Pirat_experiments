setwd("~/PhD_project/2020-proteomics-transcriptomics")

source("../data_loader/maxquant_data.R")
source("../data_loader/proline_data.R")
source("../data_loader/huang2020_data.R")
# library(reticulate)
# library(missForest)
# reticulate::use_condaenv("classics")
# reticulate::source_python("llk_maximize.py")

put_imputed_in_or = function(imputed.ds, or.ds) {
  imputed.names = colnames(imputed.ds)
  or.names = colnames(or.ds)
  is.imputed = or.names %in% imputed.names
  or.ds[, is.imputed] = imputed.ds
  return(or.ds)
}

# Load Maxquant data
# data.name = "MaxQuant"
# npcs = 2
# pathpeptides = "../../datasets/maxquant/peptides.txt"
# pathsupmat2 = "../../datasets/maxquant/mcp.M113.031591-2.csv"
# mq.data.comp = get.maxquant.data(pathpeptides, pathsupmat2, n.na.max = 7, rm.nested.pg = T)
# saveRDS(mq.data.comp, file = file.path("..", "processed_data", "cox_comp.rds"))
# mq.data.comp = readRDS(file.path("..", "processed_data", "cox_comp.rds"))
# groups = factor(rep(1:2, each = 4))

# Load Proline data
data.name = "Proline"
npcs = 10
pathpeptides = "../../datasets/Proline/proline_pep_prot.csv"
mq.data.comp = get.proline.data(pathpeptides, n.na.max = 39, rm.nested.pg = T)
saveRDS(mq.data.comp, file = file.path("..", "processed_data", "bouyssie_comp.rds"))
mq.data.comp = readRDS(file.path("..", "processed_data", "bouyssie_comp.rds"))
groups = factor(rep(1:10, each = 4))

# n.pg = ncol(mq.data.comp$adj)
# idxs.pg.rm = sample(n.pg, n.pg - 5, F)
# mq.data.comp.reduced = rm_pg_from_idx_merge_pg(mq.data.comp, idxs.pg.rm)
# mq.data.comp.reduced$peptides_ab = mq.data.comp.reduced$peptides_ab[, rowSums(mq.data.comp.reduced$adj) >= 1]
# mq.data.comp.reduced$adj = mq.data.comp.reduced$adj[rowSums(mq.data.comp.reduced$adj) >= 1, ]
# saveRDS(mq.data.comp.reduced, file.path("..", "processed_data", "proline_small_5pg.rds"))

# Load Huang2020 data
# data.name = "Huang2020"
# npcs = 5
# pathpeptides = "../../datasets/Huang2020/Spike-in-biol-var-OT-SN-Report.txt"
# # mq.data = get.huang.data(pathpeptides, n.na.max = 24, rm.nested.pg = T)
# # saveRDS(mq.data.comp, file = file.path("..", "processed_data", "Huang2020_data_comp.rds"))
# mq.data.comp = readRDS(file.path("..", "processed_data", "Huang2020_data_comp.rds"))
# groups = factor(rep(1:5, each = 5))


npeps = ncol(mq.data.comp$peptides_ab)

seednum = 543210
path2saveRDS = file.path("..", "experiments", data.name, paste(seednum, "_w_shared", sep = ""))
saveRDS(mq.data.comp, file = file.path(path2saveRDS, "DATA.rds"))
dir.create(path2saveRDS)
set.seed(seednum)


# Our method
source("pipeline_impute.R")
start_time <- Sys.time()
res.mle_mnar = pipeline_llkimpute(mq.data.comp, pep.ab.comp=NULL, nu_factor=2)
end_time = Sys.time()
saveRDS(res.mle_mnar, file = file.path(path2saveRDS, "MLEMNAR.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(path2saveRDS, "MLEMNAR_time.rds"))

# # Our method MCAR version
# source("pipeline_impute.R")
# start_time <- Sys.time()
# res.mle_mcar = pipeline_llkimpute(mq.data.comp, nu_factor=2, mcar = T)
# end_time = Sys.time()
# saveRDS(res.mle_mcar, file = file.path(path2saveRDS, "MLEMCAR.rds"))
# difference <- difftime(end_time, start_time, units='mins')
# saveRDS(difference, file.path(path2saveRDS, "MLEMCAR_time.rds"))
# 
# # Our method MCAR, transpose version
# source("pipeline_impute.R")
# start_time <- Sys.time()
# res.mle_mcar = pipeline_llkimpute(mq.data.comp, nu_factor=2, mcar = T, transpose = T)
# end_time = Sys.time()
# saveRDS(res.mle_mcar, file = file.path(path2saveRDS, "MLEMCAR_transpose.rds"))
# difference <- difftime(end_time, start_time, units='mins')
# saveRDS(difference, file.path(path2saveRDS, "MLEMCAR_transpose_time.rds"))
# 
# MsImpute MNAR
library(msImpute)
start_time <- Sys.time()
res.msimpute_mnar = t(msImpute(
  t(
    mq.data.comp$peptides_ab[, colSums(!is.na(mq.data.comp$peptides_ab)) >= 4]
  ), method = "v2-mnar", group = groups))
res.msimpute_mnar = put_imputed_in_or(res.msimpute_mnar, mq.data.comp$peptides_ab)
end_time = Sys.time()
saveRDS(res.msimpute_mnar, file = file.path(path2saveRDS, "msImpute_mnar.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(path2saveRDS, "msImpute_mnar_time.rds"))

# MLE classics
library("MsCoreUtils")
start_time <- Sys.time()
res.mleclass = impute_mle(mq.data.comp$peptides_ab)
end_time = Sys.time()
saveRDS(res.mleclass, file = file.path(path2saveRDS, "MLECLASS.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(path2saveRDS, "MLECLASS_time.rds"))


# MsImpute
library(msImpute)
start_time <- Sys.time()
res.msimpute = t(msImpute(
  t(
    mq.data.comp$peptides_ab[, colSums(!is.na(mq.data.comp$peptides_ab)) >= 4]
  ), method = "v2", group = groups))
res.msimpute = put_imputed_in_or(res.msimpute, mq.data.comp$peptides_ab)
end_time = Sys.time()
saveRDS(res.msimpute, file = file.path(path2saveRDS, "msImpute.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(path2saveRDS, "msImpute_time.rds"))
# 
# # MissForest
# library("missForest")
# start_time <- Sys.time()
# res.misfor = missForest(mq.data.comp$peptides_ab)
# end_time = Sys.time()
# saveRDS(res.misfor$ximp, file = file.path(path2saveRDS, "MissForest.rds"))
# difference <- difftime(end_time, start_time, units='mins')
# saveRDS(difference, file.path(path2saveRDS, "MissForest_time.rds"))
# 
# # KNN
# library("impute")
# start_time <- Sys.time()
# res.knn = impute.knn(t(mq.data.comp$peptides_ab), k = 10, rowmax = 1, colmax = 1)
# end_time = Sys.time()
# res.knn = t(res.knn$data)
# saveRDS(res.knn, file = file.path(path2saveRDS, "KNN.rds"))
# difference <- difftime(end_time, start_time, units='mins')
# saveRDS(difference, file.path(path2saveRDS, "KNN_time.rds"))
# 
# 
# # QRILC
# library("imputeLCMD")
# start_time <- Sys.time()
# res.qrilc = impute.QRILC(t(mq.data.comp$peptides_ab))
# end_time = Sys.time()
# res.qrilc = t(res.qrilc[[1]])
# saveRDS(res.qrilc, file = file.path(path2saveRDS, "QRILC.rds"))
# difference <- difftime(end_time, start_time, units='mins')
# saveRDS(difference, file.path(path2saveRDS, "QRILC_time.rds"))
# 
# 
# # MinProb
# library("imputeLCMD")
# start_time <- Sys.time()
# res.minprob = t(impute.MinProb(t(mq.data.comp$peptides_ab)))
# end_time = Sys.time()
# saveRDS(res.minprob, file = file.path(path2saveRDS, "MinProb.rds"))
# difference <- difftime(end_time, start_time, units='mins')
# saveRDS(difference, file.path(path2saveRDS, "MinProb_time.rds"))
# 
# # SVD
# library(pcaMethods)
# start_time <- Sys.time()
# res.svd = pca(t(mq.data.comp$peptides_ab), nPcs = 2, method = "svdImpute")
# end_time = Sys.time()
# res.svd = t(res.svd@completeObs)
# saveRDS(res.svd, file = file.path(path2saveRDS, "SVD.rds"))
# difference <- difftime(end_time, start_time, units='mins')
# saveRDS(difference, file.path(path2saveRDS, "SVD_time.rds"))
# 
# # LLS
# library(pcaMethods)
# start_time <- Sys.time()
# res.lls = llsImpute(mq.data.comp$peptides_ab)
# end_time = Sys.time()
# res.lls = res.lls@completeObs
# saveRDS(res.lls, file = file.path(path2saveRDS, "LLS.rds"))
# difference <- difftime(end_time, start_time, units='mins')
# saveRDS(difference, file.path(path2saveRDS, "LLS_time.rds"))
# 
# # trKNN
# source('Imput_funcs.r')
# library(tidyverse)
# sim_trKNN_wrapper <- function(data) {
#   result <- data %>% as.matrix %>% t %>% imputeKNN(., k=10, distance='truncation') %>% t
#   return(result)
# }
# start_time <- Sys.time()
# res.trknn <- sim_trKNN_wrapper(
#   mq.data.comp$peptides_ab[, !colSums(is.na(mq.data.comp$peptides_ab)) >= 3])
# res.trknn = put_imputed_in_or(res.trknn, mq.data.comp$peptides_ab)
# end_time = Sys.time()
# saveRDS(res.trknn, file = file.path(path2saveRDS, "trKNN.rds"))
# difference <- difftime(end_time, start_time, units='mins')
# saveRDS(difference, file.path(path2saveRDS, "trKNN_time.rds"))
# 
# # Seq-KNN
# library(SeqKnn)
# start_time <- Sys.time()
# res.seqknn <- t(SeqKNN(t(mq.data.comp$peptides_ab), k = 10))
# end_time = Sys.time()
# saveRDS(res.seqknn, file = file.path(path2saveRDS, "seqKNN.rds"))
# difference <- difftime(end_time, start_time, units='mins')
# saveRDS(difference, file.path(path2saveRDS, "seqKNN_time.rds"))
# 
# 
# # GMS
# library(GMSimpute)
# start_time <- Sys.time()
# res.gms <- t(GMS.Lasso(
#   t(
#     mq.data.comp$peptides_ab[, colSums(!is.na(mq.data.comp$peptides_ab)) >= 3]
#     ),nfolds=3,log.scale=FALSE,TS.Lasso=TRUE))
# res.gms = put_imputed_in_or(res.gms, mq.data.comp$peptides_ab)
# end_time = Sys.time()
# saveRDS(res.gms, file = file.path(path2saveRDS, "GMS.rds"))
# difference <- difftime(end_time, start_time, units='mins')
# saveRDS(difference, file.path(path2saveRDS, "GMS_time.rds"))
# 
# 
# Impseq
library(rrcovNA)
start_time <- Sys.time()
res.impseq <- t(impSeq(t(mq.data.comp$peptides_ab)))
end_time = Sys.time()
saveRDS(res.impseq, file = file.path(path2saveRDS, "ImpSeq.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(path2saveRDS, "ImpSeq_time.rds"))
# 
# 
# # ImpSeqRob
library(rrcovNA)
start_time <- Sys.time()
res.impseqrob <- t(impSeqRob(t(mq.data.comp$peptides_ab))$x)
end_time = Sys.time()
saveRDS(res.impseqrob, file = file.path(path2saveRDS, "ImpSeqRob.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(path2saveRDS, "ImpSeqRob_time.rds"))
# 
# 
# # BPCA
library("pcaMethods")
start_time <- Sys.time()
res.bpca = pca(t(mq.data.comp$peptides_ab), method = "bpca", nPcs = npcs)
end_time = Sys.time()
res.bpca = t(res.bpca@completeObs)
saveRDS(res.bpca, file = file.path(path2saveRDS, "BPCA.rds"))
difference <- difftime(end_time, start_time, units='mins')
saveRDS(difference, file.path(path2saveRDS, "BPCA_time.rds"))
