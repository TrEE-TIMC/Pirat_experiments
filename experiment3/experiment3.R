setwd("~/PhD_project/2020-proteomics-transcriptomics/experiment_script/")

source("../data_loader/reset_data.R")
source("../method/utils.R")
library(reticulate)
library(foreach)
library(doParallel)
reticulate::use_condaenv("classics")
reticulate::source_python("llk_maximize.py")

put_imputed_in_or = function(imputed.ds, or.ds) {
  imputed.names = colnames(imputed.ds)
  or.names = colnames(or.ds)
  is.imputed = or.names %in% imputed.names
  or.ds[, is.imputed] = imputed.ds
  return(or.ds)
}

# Load RESET data
# name.data = "RESET"
# cond_idx_rna = 1:18
# cond_idx_pep = 1:18
# npcs = 6
# path = "../datasets/reset_RNA_prot.txt"
# 
# sel_pep = c("W1_1_MS", "W1_2_MS", "W1_3_MS","W2_1_MS", "W2_2_MS", "W2_3_MS",
#             "R1_1_MS", "R1_2_MS", "R1_3_MS", "R2_1_MS", "R2_2_MS", "R2_3_MS",
#             "R3_1_MS", "R3_2_MS", "R3_3_MS", "R4_1_MS", "R4_2_MS", "R4_3_MS")
# sel_rna = c("W1_1_RNA", "W1_2_RNA", "W1_3_RNA", "W2_1_RNA", "W2_2_RNA", "W2_3_RNA",
#             "R1_1_RNA", "R1_2_RNA", "R1_3_RNA", "R2_1_RNA", "R2_2_RNA", "R2_3_RNA",
#             "R3_1_RNA", "R3_2_RNA", "R3_3_RNA", "R4_1_RNA", "R4_2_RNA", "R4_3_RNA")


# Load MouseColon data
name.data = "MouseColon"
npcs = 6
cond_idx_rna = c(1,1,1,2,2,2,2,3,3,3,3,3,4,4,4,5,5,6,6,6,6,6)
cond_idx_pep = rep(1:6, each=3)
path = "../datasets/mousecolon_RNA_prot.txt"

sel_pep = c("AbsProBR1_MS", "AbsProBR2_MS", "AbsProBR3_MS", "SecPDGBR1_MS",
            "SecPDGBR2_MS", "SecPDGBR3_MS", "TuftBR1_MS", "TuftBR2_MS",
            "TuftBR3_MS", "StemBR1_MS", "StemBR2_MS", "StemBR3_MS", "EECBR1_MS",
            "EECBR2_MS", "EECBR3_MS", "EntBR1_MS", "EntBR2_MS", "EntBR3_MS")

sel_rna = c("AbsProBR1_RNAseq", "AbsProBR2_RNAseq", "AbsProBR3_RNAseq",
            "SecPDGBR1_RNAseq", "SecPDGBR2_RNAseq", "SecPDGBR3_RNAseq", "SecPDGBR4_RNAseq",
            "TuftBR1_RNAseq", "TuftBR2_RNAseq", "TuftBR3_RNAseq", "TuftBR4_RNAseq", "TuftBR5_RNAseq",
            "StemBR1_RNAseq", "StemBR2_RNAseq", "StemBR3_RNAseq",
            "EECBR1_RNAseq", "EECBR2_RNAseq", "EnterocyteBR1_RNAseq",
            "EnterocyteBR2_RNAseq", "EnterocyteBR3_RNAseq", "EnterocyteBR4_RNAseq",
            "EnterocyteBR5_RNAseq")

# all.data = get.reset.data(path, sel_pep, sel_rna)
# 
# nsamples = nrow(all.data$peptides_ab)
# cat("Size raw dataset", dim(all.data$adj))
# 
# all.data.wo.bad.rna = remove_Os_rna_reset(all.data, cond_idx_rna)
# dim(all.data.wo.bad.rna$rnas_ab)
# table(colSums(all.data.wo.bad.rna$rnas_ab == 0))
# 
# tol_na_pep = nsamples - 1  # Remove peptides that have 'tol_na_pep' or more Missing values
# all.data.comp = remove_NA_pep_reset(all.data.wo.bad.rna, tol_na_pep)
# 
# # Remove nested prots x1
# print("Remove nested prots x1...")
# idx.emb.prots = get_indexes_embedded_prots(all.data.comp$adj)
# all.data.comp = rm_pg_from_idx_merge_pg(all.data.comp, idx.emb.prots)
# 
# saveRDS(all.data.comp, "all_data_comp_pg_RESET.rds")
all.data.comp = readRDS("../processed_data/all_data_comp_pg_MouseColon.rds")


nsamples = nrow(all.data.comp$peptides_ab)

npeps = ncol(all.data.comp$peptides_ab)

prop.na = mean(is.na(all.data.comp$peptides_ab))
cat("NA proportion", prop.na)
cat("N shared pap", sum(rowSums(all.data.comp$adj) >= 2))
cat("Size preprocessed dataset", dim(all.data.comp$adj))


# Add pseudo missing values
tol.pseudo.na = nsamples - 1
n.na = prop.na * npeps * nsamples
n.pseudo.na = 0.01 * n.na
pseudo.mv.rate = n.pseudo.na / ((1 - prop.na) * npeps * nsamples)
seednums = c(543213:543214) #  c(543210, 543211, 543212, 543213, 543214, 543215, 543216, 543217, 543218, 543219)
mnar.mv.rates = c(1)
mncar.fold = paste("MNAR", mnar.mv.rates)
sd.t = sd(all.data.comp$peptides_ab, na.rm = T)/2
# Get approximated function of censoring mechanism mean vs given MNAR proportion
expec_sum_m = Vectorize(
  function(mu, xx, sd) {return(1 - mean(pnorm(xx, mu, sd), na.rm = T))},
  vectorize.args = "mu")
qlow = quantile(all.data.comp$peptides_ab, probs = pseudo.mv.rate*0.001, na.rm = T)
qhigh = quantile(all.data.comp$peptides_ab, probs = pseudo.mv.rate, na.rm = T)
q_vec = seq(qlow, qhigh, length.out = 100)
expec_sum_m_vec = expec_sum_m(q_vec, all.data.comp$peptides_ab, sd.t)
approx_q = approxfun(expec_sum_m_vec, q_vec)

# n.cores = 2
# cur.cluster = makeCluster(n.cores, type = "FORK", outfile = "")
# print(cur.cluster)
# registerDoParallel(cl = cur.cluster)
# getDoParRegistered()
# getDoParWorkers()

registerDoSEQ()

n.seeds = length(seednums)
n.mnars = length(mnar.mv.rates)

impmethods = c("DATA", "ImpSeq", "KNN", "BPCA", "MLECLASS", "GMS",
               "QRILC", "MinProb", "SVD", "LLS", "trKNN", "seqKNN", "ImpSeqRob") # Imputation methods

impmethods = c("MLEMNAR_PG")

overall_start <- Sys.time()
foreach(seednum = seednums) %:%
  foreach(idx.mnar.mv.rate = 1:n.mnars) %:%
    foreach(method = impmethods) %dopar% {
      print("###################")
      print(paste(seednum, idx.mnar.mv.rate, method))
      print("###################")
      mnar.mv.rate = mnar.mv.rates[idx.mnar.mv.rate]
      path2saveRDS = file.path("../experiments", name.data, 
                               "pseudo_na_percent_1_mis", 
                               mncar.fold[idx.mnar.mv.rate])
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
      
      # Remove peptide w too many NAs from pseudo dataset in both control and pseudo datasets
      idxpep.fullna = which(colSums(is.na(all.pseudo$peptides_ab)) >= tol.pseudo.na)
      reset.comp.rmfullna = remove_pep_from_idx(all.data.comp, idxpep.fullna)
      all.pseudo = remove_NA_pep_reset(all.pseudo, tol.pseudo.na)
      all.pseudo[["comp_pep_abs"]] = reset.comp.rmfullna$peptides_ab
      
      print(paste("NA proportion", mean(is.na(all.pseudo$peptides_ab))))
      print(paste("Dim adj", paste(dim(all.pseudo$adj))))
      print(paste("N shared pep", sum(rowSums(all.pseudo$adj) >= 2)))
      print(paste("N pseudo NA", sum(is.na(all.pseudo$peptides_ab) & !is.na(all.pseudo$comp_pep_abs))))
      
      if (method == "DATA") {
        print("Remove nested prots x2...")
        idx.emb.prots = get_indexes_embedded_prots(all.pseudo$adj)
        all.pseudo = rm_pg_from_idx_merge_pg(all.pseudo, idx.emb.prots)
        saveRDS(all.pseudo, file = file.path(path2saveRDS, "DATA.rds"))
      }
      
      
      # Our method
      if (method == "MLEMNAR") {
        source("../method/pipeline_impute.R")
        start_time <- Sys.time()
        res.mle_mnar = pipeline_llkimpute(all.pseudo, 
                                          # pep.ab.comp = all.pseudo$comp_pep_abs,
                                          nu_factor = 2)
        end_time = Sys.time()
        saveRDS(res.mle_mnar, file = file.path(path2saveRDS, "MLEMNAR_full.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "MLEMNAR_full_time.rds"))
      }
      
      # Our method Transposed for PG of size 1
      if (method == "MLEMNAR_TPG1") {
        source("../method/pipeline_impute.R")
        start_time <- Sys.time()
        path.normal.imputed = file.path(path2saveRDS, "MLEMNAR_full.rds") 
        res.mle_mnar.tpg1 = pipeline_llkimpute(all.pseudo, 
                                               pep.ab.comp = all.pseudo$comp_pep_abs,
                                               nu_factor = 2,
                                               pathifcc1 = path.normal.imputed)
        end_time = Sys.time()
        saveRDS(res.mle_mnar.tpg1, file = file.path(path2saveRDS, "MLEMNAR_TPG1.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "MLEMNAR_TPG1_time.rds"))
      }
      
      # Our method ProteoGenomics
      if (method == "MLEMNAR_PG") {
        source("../method/pipeline_impute.R")
        start_time <- Sys.time()
        res.mle_mnar_pg = pipeline_llkimpute(all.pseudo,
                                          pep.ab.comp = all.pseudo$comp_pep_abs,
                                          nu_factor = 2,
                                          rna.cond.mask = cond_idx_rna,
                                          pep.cond.mask = cond_idx_pep,
                                          protidxs = NULL,
                                          max.pg.size2imp = 1)
        end_time = Sys.time()
        saveRDS(res.mle_mnar_pg, file = file.path(path2saveRDS, "MLEMNAR_PG.rds"))
        difference <- difftime(end_time, start_time, units='mins')
        saveRDS(difference, file.path(path2saveRDS, "MLEMNAR_PG_time.rds"))
      }
      
      print(paste(method, "DONE !"))
      NULL
    }
overall_end <- Sys.time()
stopCluster(cl = cur.cluster)
print(overall_end - overall_start)
