source("PEMM_fun.R")
library(Pirat)
library(ggplot2)
data("ropers")

#############################
### Test PEMM on Ropers2021
#############################

set.seed(10)
obs2NApep <- ropers$peptides_ab[,colSums(is.na(ropers$peptides_ab)) <= 0]
psi.df = estimate_psi_df(obs2NApep)
gammas = estimate_gamma(ropers$peptides_ab)
has_converged = rep(T, ncol(ropers$adj))
var_exploded = rep(F, ncol(ropers$adj))
maxiter_reached = rep(F, ncol(ropers$adj))
llk_is_na = rep(F, ncol(ropers$adj))
llk_decreased = rep(F, ncol(ropers$adj))
for (idx.pg in 1:ncol(ropers$adj)) {#ncol(ropers$adj)) {
  print(idx.pg)
  pep.abs = ropers$peptides_ab[, ropers$adj[, idx.pg] == T]
  if (sum(ropers$adj[, idx.pg]) > 1) {
    if (all(rowSums(!is.na(pep.abs)) != 0)) {
      lambda = psi.df$psi
      K = ncol(pep.abs) + 2*psi.df$df
      resu = PEMM_fun(pep.abs, gammas$gamma_1, lambda = lambda, K = K, pos = T,
                      tol = 0.0001,
                     maxIter = 10000)
      if (any(resu$Sigma >= 1000)) {
        print("Big VAR")
        var_exploded[idx.pg] = T
      }
      if (!resu$has_converged) {
        print("Doesn't converge")
        has_converged[idx.pg] = resu$has_converged
      }
      if (is.unsorted(na.omit(resu$illik[-1]))) {
        print("Pen LLK decreases")
        llk_decreased[idx.pg] = T
      }
      maxiter_reached[idx.pg] = resu$mxiter_reached
      llk_is_na[idx.pg] = resu$llk_is_na
    }
  }
}
sum(has_converged)
sum(var_exploded)
sum(maxiter_reached)
sum(llk_is_na)
sum(llk_decreased)

# PG 34
idx.pg = 34
pep.abs = ropers$peptides_ab[, ropers$adj[, idx.pg] == T]
K = ncol(pep.abs) + 2*psi.df$df
resu = PEMM_fun(pep.abs, gammas$gamma_1, lambda = lambda, K = K, pos = T,
                tol = 0.0001, maxIter = 10000)
resu$illiks
diag(resu$Sigma)
plot(resu$illiks)


###########
# Bias on mu with synthetic data
###########

reticulate::use_condaenv("r-pirat")
py <- reticulate::import("PyPirat", delay_load = TRUE)

n.datasets = 1000
rmses.Pirat = c()
rmses.PEMM = c()
bias.Pirat = c()
bias.PEMM = c()
p = 10
n = 10
set.seed(1234)
s1 = 25
s2 = 35
rho = 0.5
gamma1 = 0.5
sigma = toeplitz(rho^(0:(p-1)))
has_converged = rep(T, n.datasets)
var_exploded = rep(F, n.datasets)
maxiter_reached = rep(F, n.datasets)
llk_is_na = rep(F, n.datasets)
mu_exploded = rep(F, n.datasets)
i = 1
while(i <= n.datasets) {
  print(i)
  mu = runif(p, s1, s2)
  X = matrix(rnorm(n * p), n) %*% chol(sigma)
  X = t(t(X)+mu)
  gamma0 = - gamma1 * s1
  p.mis = pmin(1, exp(-gamma0 - gamma1*X))
  M = matrix(rbinom(n*p, 1, c(p.mis)), n)
  M = M == 1
  X.w.na = X
  X.w.na[M] = NA
  if (all(rowSums(!is.na(X.w.na)) > 0) & all(colSums(!is.na(X.w.na)) > 1)) {
    alpha = 1.
    beta = 1.
    resu.PEMM = tryCatch(PEMM_fun(X.w.na, gamma1, lambda = beta, K = alpha + p,
                                  pos = T, tol = 0.0001, maxIter = 1000), 
                         error = function(e) {message("Error occurred: ", e); return(list())})
    
    if (length(resu.PEMM) != 0) {
      if (any(resu.PEMM$Sigma >= 1000)) {
        var_exploded[i] = T
      }
      if (any(abs(resu.PEMM$mu) >= 1000)) {
        mu_exploded[i] = T
      }
      if (!resu.PEMM$has_converged) {
        has_converged[i] = resu.PEMM$has_converged
      }
      maxiter_reached[i] = resu.PEMM$mxiter_reached
      llk_is_na[i] = resu.PEMM$llk_is_na
      
      mu.PEMM = resu.PEMM$mu
      bias.PEMM = c(bias.PEMM, mean(mu.PEMM - mu))  
      rmse.mu.PEMM = sqrt(mean((mu.PEMM - mu)^2))
      rmses.PEMM = c(rmses.PEMM, rmse.mu.PEMM)

      res.Pirat = py$estimate_params_and_impute(X.w.na, phi0 = gamma0,
                                                phi = gamma1,
                                                K = (alpha + p) * 2,
                                                psi = beta * 2,
                                                eps_chol = 1e-4,
                                                eps_phi = 1e-5,
                                                tol_obj = 1e-7,
                                                tol_grad = 1e-5,
                                                tol_param = 1e-4,
                                                maxiter = as.integer(10000),
                                                lr = 0.5,
                                                phi_known = TRUE,
                                                max_try = 50,
                                                max_ls = 500,
                                                eps_sig = 1e-4,
                                                nsamples = 1000)
      mu.Pirat = res.Pirat$mu
      rmse.mu.Pirat = sqrt(mean((mu.Pirat - mu)^2))
      rmses.Pirat = c(rmses.Pirat, rmse.mu.Pirat)
      bias.Pirat = c(bias.Pirat, mean(mu.Pirat - mu))
      i = i + 1
    }
  }
}

sum(!has_converged)
sum(maxiter_reached)
sum(llk_is_na)
sum(var_exploded)
sum(mu_exploded)

data.vis = data.frame(rmses = c(rmses.PEMM, rmses.Pirat),
                      bias = c(bias.PEMM, bias.Pirat),
                      has_converged = c(has_converged, has_converged),
                      method = factor(c(rep(c("PEMM", "Pirat"), each = n.datasets))))
saveRDS(data.vis, "data.vis")
data.vis1 = data.vis[data.vis$has_converged, ]
g = ggplot(data.vis1, aes(x = method, y = rmses)) + geom_boxplot() # + ylim(c(0.1, 10000))
g <- g + theme(panel.background = element_blank()) + ylab("RMSE") + xlab("Method")
g
ggsave("rmses_gamma1_05.pdf",
       width = 4, height = 4,
       device = "pdf")

g = ggplot(data.vis1, aes(x = method, y = bias)) + geom_boxplot() # + ylim(c(0.1, 10000))
g <- g + theme(panel.background = element_blank()) + ylab("Mean bias") + xlab("Method")
g
ggsave("bias_gamma1_05.pdf",
       width = 4, height = 4,
       device = "pdf")
wilcox.test(data.vis1$bias[data.vis1$method == "PEMM"], data.vis1$bias[data.vis1$method == "Pirat"], alternative = c("two.sided"))
wilcox.test(data.vis1$rmses[data.vis1$method == "Pirat"], data.vis1$rmses[data.vis1$method == "PEMM"], alternative = c("two.sided"))
