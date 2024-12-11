# llk_maximize.py

import numpy as np
import torch as t
#from PyPirat import FullBatchLBFGS
#from LBFGS import FullBatchLBFGS
import time
import PyPirat




def get_sigma_from_log_chol(log_chol_m, eps_chol=1e-10, eps_sig=1e-4):
    eps_chol = t.tensor(eps_chol)
    p = log_chol_m.size(0)
    diag_chol = t.diag(log_chol_m)
    log_chol_copy = log_chol_m.clone()
    log_chol_copy[t.arange(p), t.arange(p)] = t.exp(diag_chol) + eps_chol
    # print(diag_chol)
    # print(t.diag(log_chol_copy))
    # print(t.sum(t.linalg.inv(log_chol_copy @ log_chol_copy.t())))
    return log_chol_copy @ log_chol_copy.t() + eps_sig*t.eye(p)


def compute_llk(X, mu, sigma, phi0=0, phi=0, nsamples=10000, method="low_bound"):
    llik = 0
    # prob_max = t.tensor(1 - 1e-6)
    npep = X.size(1)
    M = t.isnan(X)
    uniq_m_patterns, idx_uniq = t.unique(M, dim=0, return_inverse=True)
    n_uniq = uniq_m_patterns.size(0)
    for i_uniq in range(n_uniq):
        row_ret = t.where(idx_uniq == i_uniq)[0]
        Xi = X[row_ret, :]
        idxi = t.where(~uniq_m_patterns[i_uniq])[0]
        Xi = Xi[:, idxi]
        Si = sigma[idxi, :][:, idxi]
        nmis = npep - idxi.size(0)
        try:
            Sii_Xi_tilde = t.linalg.solve(Si, (Xi - mu[idxi]).t())
        except RuntimeError:
            Sii_Xi_tilde = t.nan * t.ones((Xi - mu[idxi]).t().size(), dtype=t.float64)
            nmis = 0
        log_det_Si = t.logdet(Si)
        oo = -row_ret.size(0)*log_det_Si - t.sum((Xi - mu[idxi]) * Sii_Xi_tilde.t())
        oo = 0.5 * oo
        if phi == 0.:
            vv = 0
        else:
            # vv1 = t.sum(t.log(1 - t.min(t.exp(-phi0 - phi*Xi), prob_max)))
            if nmis == 0:
                vv3_exact = 0
            else:
                idxmi = t.where(uniq_m_patterns[i_uniq])[0]
                Smi = sigma[idxmi, :][:, idxmi]
                Smi_xi = sigma[idxmi, :][:, idxi]
                mu_mis = mu[idxmi] + t.matmul(Smi_xi, Sii_Xi_tilde).t() # nret x nmis
                S_mis = Smi - t.matmul(Smi_xi, t.linalg.solve(Si, Smi_xi.t()))

                ## Exact method 1
                if method == "approx":
                    sigmas_mi = t.diag(S_mis)
                    normal = t.distributions.normal.Normal(0, 1)
                    if t.any(sigmas_mi <= 0) or t.any(t.isnan(sigmas_mi)) or t.any(t.isinf(sigmas_mi)):
                        vv3_exact = t.nan
                    else:
                        vv3_exact = normal.cdf((-phi0/phi - mu_mis)/t.sqrt(sigmas_mi)) +\
                                    t.exp(-phi0 - phi*mu_mis + 0.5*phi**2*sigmas_mi) *\
                                    (1 - normal.cdf((-phi0/phi - mu_mis + phi*sigmas_mi)/t.sqrt(sigmas_mi)))
                        vv3_exact = t.sum(t.log(vv3_exact))

                ##Exact method 2
                elif method == "low_bound":
                    sigmas_mi = t.diag(S_mis)
                    normal = t.distributions.normal.Normal(0, 1)
                    if t.any(sigmas_mi <= 0) or t.any(t.isnan(sigmas_mi)) or t.any(t.isinf(sigmas_mi)) or \
                            t.any(t.isnan(mu_mis)):
                        vv3_exact = t.nan
                    else:
                        sd_sig_mi = t.sqrt(sigmas_mi)
                        vv3_exact = -(1 - normal.cdf((-phi0/phi - mu_mis)/sd_sig_mi)) * (phi0 + phi*mu_mis) - \
                            phi * sd_sig_mi * t.exp(normal.log_prob((-phi0/phi - mu_mis)/sd_sig_mi))
                        vv3_exact = t.sum(vv3_exact)

                ## Sample method
                elif method == "sampling":
                    samples = t.normal(0, 1, size=(row_ret.size(0), nmis, int(nsamples)), dtype=t.float64) # nret x nmis x nsamples
                    # try:
                    #     chol_S_mis = t.linalg.cholesky(S_mis)
                    # except RuntimeError:
                    #     chol_S_mis = t.nan * t.ones(S_mis.size(), dtype=t.float64)
                    chol_S_mis = t.linalg.cholesky_ex(S_mis)[0]
                    Xmis_samples = chol_S_mis @ samples +\
                                   mu_mis.unsqueeze(2).repeat(1, 1, int(nsamples))  # nret x nmis x nsamples
                    p_mis = t.clamp(t.exp(-phi0 - phi * Xmis_samples), max=1.)
                    vv3_exact = t.sum(t.log(t.mean(t.prod(p_mis, dim=1), 1)))
            vv = vv3_exact
        llik += oo + vv
    return llik

def compute_pen_llk(X, mu, log_chol_params, K, psi, phi0, phi, eps_chol=1e-5, eps_sig=1e-6, nsamples=1000,
                    method="low_bound"):
    p = X.size(1)
    # mask = t.tril(t.ones(p, p)) # Enable to autograd only lower tri matrix
    idx_ltri = t.tril_indices(p, p)
    log_chol_m = t.zeros((p, p), dtype=t.float64)
    log_chol_m[idx_ltri[0], idx_ltri[1]] = log_chol_params
    sigma = get_sigma_from_log_chol(log_chol_m, eps_chol=eps_chol, eps_sig=eps_sig)
    try:
        pen_trace = - 0.5 * t.trace(t.linalg.solve(sigma.t(), psi.t()))
    except RuntimeError:
        pen_trace = t.nan
    llk = compute_llk(X, mu, sigma, phi0, phi, nsamples, method=method)
    pen_det = - K / 2. * t.logdet(sigma)
    return llk + pen_det + pen_trace


def get_matrix_from_tri_coefs(sigma_tri, input_low=True, out_low=True, out_up=False):
    p = int((-1 + np.sqrt(1. + 8*sigma_tri.shape[0]))/2)
    sigma = np.zeros([p, p])
    if input_low:
        uidx = np.tril_indices(p)
        if out_low:
            sigma[uidx] = sigma_tri
        if out_up:
            sigma.transpose()[uidx] = sigma_tri
    else:
        uidx = np.triu_indices(p)
        if out_up:
            sigma[uidx] = sigma_tri
        if out_low:
            sigma.transpose()[uidx] = sigma_tri
    return sigma


def impute_from_params(X_w_na, mu, sigma, phi, phi0):
    X = np.copy(X_w_na)
    npep = X.shape[1]
    var_imp = np.zeros(X.shape)
    mu = np.asarray(mu)
    M = np.isnan(X)
    nsample = 1000
    uniq_m_patterns, idx_uniq = np.unique(M, axis=0, return_inverse=True)
    n_uniq = uniq_m_patterns.shape[0]
    for i_uniq in range(n_uniq):
        ii = np.where(uniq_m_patterns[i_uniq])[0]
        if ii.shape[0] >= 1:
            row_ret = np.where(idx_uniq == i_uniq)[0]
            ii_obs = np.where(~uniq_m_patterns[i_uniq])[0]
            Soo = sigma[np.ix_(ii_obs, ii_obs)]
            Soo_inv = np.linalg.inv(Soo)  # nobs x nobs
            mu_mis = mu[ii] + np.matmul(
                sigma[np.ix_(ii, ii_obs)], np.matmul(
                    Soo_inv, np.transpose(
                        X[np.ix_(row_ret, ii_obs)] - mu[ii_obs]))).transpose()  # size_uniq x nmis
            cov_mis = sigma[np.ix_(ii, ii)] - np.matmul(sigma[np.ix_(ii, ii_obs)],
                                                    np.matmul(Soo_inv, sigma[np.ix_(ii_obs, ii)]))  # nmis x nmis
            if phi != 0:
                nmis = npep - ii_obs.shape[0]
                ## Exact method
                # sigmas_mi = np.diag(cov_mis)
                # weird_exp = np.exp(-phi0 - phi*mu_mis + 0.5*phi**2*sigmas_mi)
                # mu_mis_tilde = mu_mis - phi*sigmas_mi
                # beta_l = (-phi0/phi - mu_mis) / sigmas_mi
                # alpha_h = (-phi0/phi - mu_mis_tilde) / sigmas_mi
                # ZZ_l = norm.cdf(beta_l)
                # ZZ_h = 1 - norm.cdf(alpha_h)
                # moment0 = norm.cdf(-phi0/phi, mu_mis, np.sqrt(sigmas_mi)) +\
                #     weird_exp * (1 - norm.cdf(-phi0/phi, mu_mis_tilde, np.sqrt(sigmas_mi)))
                # # ratio_pdfbeta_Z_l = norm.pdf(beta_l) / ZZ_l
                # # ratio_pdfbeta_Z_l[np.isnan(ratio_pdfbeta_Z_l)] = - beta_l - 1/beta_l + 2/beta_l**3  # Taylor approximation
                # # moment1 = mu_mis - norm.pdf(beta_l) / ZZ_l * sigmas_mi + weird_exp * \
                # #     (mu_mis_tilde + norm.pdf(alpha_h) / ZZ_h)
                # # moment2_l = mu_mis**2 - 2*mu_mis*sigmas_mi*norm.pdf(beta_l)/ZZ_l + \
                # #     sigmas_mi**2 * (1 - beta_l * norm.pdf(beta_l) / ZZ_l)
                # # moment2_h = mu_mis_tilde ** 2 + 2 * mu_mis_tilde * sigmas_mi * norm.pdf(alpha_h) / ZZ_h + \
                # #     sigmas_mi ** 2 * (alpha_h * norm.pdf(alpha_h) / ZZ_h + 1)
                # # moment2 = moment2_l + weird_exp*moment2_h
                # clip_a_low = -np.inf
                # clip_b_high = np.inf
                # n_row, n_col = mu_mis.shape[0], mu_mis.shape[1]
                # moment1, moment2 = np.zeros((n_row, n_col)), np.zeros((n_row, n_col))
                # for i in range(n_row):
                #     for j in range(n_col):
                #         mu_mis_ij = mu_mis[i, j]
                #         mu_mis_tilde_ij = mu_mis_tilde[i, j]
                #         sigmas_mi_j = sigmas_mi[j]
                #         weird_exp_ij = weird_exp[i, j]
                #         clib_b_low = (-phi0 / phi - mu_mis_ij) / np.sqrt(sigmas_mi_j)
                #         clip_a_high = (-phi0 / phi - mu_mis_tilde_ij) / np.sqrt(sigmas_mi_j)
                #         norm_l = norm.cdf(-phi0 / phi, mu_mis_ij, np.sqrt(sigmas_mi_j))
                #         norm_h = 1 - norm.cdf(-phi0 / phi, mu_mis_tilde_ij, np.sqrt(sigmas_mi_j))
                #         moment1[i, j] = norm_l * \
                #             truncnorm.moment(1, clip_a_low, clib_b_low, mu_mis_ij, np.sqrt(sigmas_mi_j)) +\
                #             weird_exp_ij * norm_h * truncnorm.moment(1, clip_a_high,
                #                                                      clip_b_high, mu_mis_tilde_ij, np.sqrt(sigmas_mi_j))
                #         moment2[i, j] = norm_l * \
                #             truncnorm.moment(2, clip_a_low, clib_b_low, mu_mis_ij, np.sqrt(sigmas_mi_j)) + \
                #             weird_exp_ij * norm_h * truncnorm.moment(2, clip_a_high,
                #                                                      clip_b_high, mu_mis_tilde_ij, np.sqrt(sigmas_mi_j))
                #
                # X[np.ix_(row_ret, ii)] = moment1 / moment0
                # var_imp[np.ix_(row_ret, ii)] = moment2 / moment0 - (moment1 / moment0)**2


                ## Sample method
                samples = np.random.normal(0, 1, size=(row_ret.shape[0], nmis, nsample))  # size_uniq x nmis x nsamples
                Xmis_samples = np.matmul(np.sqrt(np.array(np.diag(cov_mis), ndmin=2)), samples) + \
                    np.repeat(np.expand_dims(mu_mis, axis=2), repeats=nsample, axis=2)  # size_uniq x nmis x nsamples
                p_mis = np.minimum(np.exp(-phi0 - phi * Xmis_samples), 1.)
                moment0 = np.mean(p_mis, axis=2)  # size_uniq x nmis
                moment1 = np.mean(Xmis_samples * p_mis, 2)  # size_uniq x nmis
                moment2 =\
                    np.mean(Xmis_samples ** 2 * p_mis, 2)  # size_uniq x nmis
                mean_Ximp = moment1 / moment0
                var_Ximp = (moment2 / moment0) - mean_Ximp**2
                var_imp[np.ix_(row_ret, ii)] = var_Ximp
                # X[np.ix_(row_ret, ii)] = np.random.normal(mean_Ximp, scale=np.sqrt(var_Ximp))
                X[np.ix_(row_ret, ii)] = mean_Ximp
            else:
                var_imp[np.ix_(row_ret, ii)] = np.diag(cov_mis)
                # X[np.ix_(row_ret, ii)] = np.random.normal(mu_mis, scale=np.sqrt(var_Ximp))
                X[np.ix_(row_ret, ii)] = mu_mis
    return X, var_imp


def estimate_params_and_impute(X, phi0=None, phi=None, K=5, psi=1., phi_known=True, eps_chol=1e-4, eps_phi=1e-8,
                               tol_obj=1e-9, tol_grad=1e-5, tol_param=1e-6, maxiter=500, lr=1., true_mu=None,
                               true_sigma=None, true_X=None, verbose=False, max_try=10, max_ls=10, eps_sig=1e-5,
                               nsamples=1000):
    """
    Estimates feature mean and covariance matrix of a dataset given missingness parameters, 
    and imputes missing values accordingly.
    
    Inputs:
        X (numpy array): dataset to impute, with features in columns and samples in rows
        phi0 (float): intercept of the missingness parameters (default: None)
        phi (float): slope of the missingness parameters (default: None)
        K (float): degree of freedom in the inverse-Wishart penalty (default: 5)
        psi (float): scale factor in the inverse-Wishart penalty (default: 1)
        phi_known (bool): are missingness parameters known (default: True)
        eps_chol (float): amount added to Cholesky's factorization's diagonal to avoid ill conditionned matrix (default: 1e-4)
        eps_phi (float): amount added to phi in case it is estimated along gaussian parameters, such that phi remains strictly positive (default: 1e-8)
        tol_obj (float): tolerance on objective function during optimization (default: 1e-9)
        tol_grad (float): tolerance on gradient norm during optimization function (default: 1e-5)
        tol_param (float): tolerance on parameters norm during optimization function (default: 1e-6)
        maxiter (int): maximum number of iterations during optimization (default: 500)
        lr (float): initial learning rate (or step size) of L-BFGS (default: 1.)
        true_mu (numpy array): true value of mu (default: None)
        true_sigma (numpy array): true value of sigma (default: None)
        true_X (numpy array): true value of X (default: None)
        verbose (bool): should verbose mode be activated (default: False)
        max_try (int): maximum number of attempts if numerical instabilities appear in optmization due to determinant approximation (default: 10)
        max_ls (int): maximum number of line search attempts in L-BFGS (default: 10)
        eps_sig (float): amount added to sigma's diagonal to avoid ill conditionned matrix (default: 1e-5)
        nsamples (int): number of samples used for MOnte-Carlo at imputation step (default: 1000)
        
    
    Outputs:
      res_dic (dict): a dictionnary containing:
            . 'Xhat' (numpy array): the imputed dataset
            . 'varXhat' (numpy array): the conditional variance of imputed values, 
            . 'phi0' (float): estimated or used phi0
            . 'phi' (float): estimated or used phi
            . 'llks' (numpy array): list of likelihood values across L-BFGS iterations
            . 'mu' (numpy array): estimated mu
            . 'sigma' (numpy array): estimated sigma
            . 'error_msg' (str): error message if error occurs
            . 'reason' (str): stop criterion of optimization process
            . 'mu_errors' (numpy array): optionnal. mu errors across iterations.
            . 'sigma_errors' (numpy array): optionnal. sigma errors across iterations.
            . 'ximp_errors' (numpy array): optionnal. imputation errors across iterations.
            . 'Xgt' (numpy array): optionnal. ground truth dataset.
    
    Notes:
      . Pirat is dedicated for single imputation, thus the seed if fixed for the Monte-Carlo method used for imputation by conditionnal mean,
    
    """
                                 
    # Added by Sam 
    TRUEFAIL = None
    
    #print('Inside estimate_params_and_impute')
    t.manual_seed(12345)
    np.random.seed(12345)

    n, p = X.shape

    if isinstance(psi, (int, float)):
        psi = psi * np.eye(p)
    psi = t.tensor(psi)

    idx_pep_mis = np.where(np.sum(np.isnan(X), 0) >= 1)[0]
    mu0 = np.nanmean(X, 0)

    # Initiate initial values of phi
    if phi is None:
        phi_init = t.tensor(0., dtype=t.float64)
    else:
        if phi_known:
            phi_init = t.tensor(phi, dtype=t.float64)
        else:
            log_phi_init = t.tensor(np.log(phi), dtype=t.float64)
            # Size must be non zero for optimizer
    if phi0 is None:
        phi0_init = t.tensor(-np.log(np.mean(1. * np.isnan(X))), dtype=t.float64)
    else:
        phi0_init = t.tensor(phi0, dtype=t.float64)

    X = t.tensor(X, dtype=t.float64)

    out = False
    lr_init = lr
    func_evals = 0
    interpolate = False
    inplace = False
    line_search = "Armijo"
    ls_debug = False
    damping = True
    eta = 2
    time_cur = time.process_time()
    # main loop
    nb_try = 0
    error_msg = "loss_eq_nan"
    while (nb_try < max_try) and (error_msg in ["var_exploded", "loss_eq_nan", "ill_posed"]):
        nb_try += 1

        # Get log-Cholesky parametrization
        # spread_unif = 0.01
        sigma0 = psi.numpy()*np.eye(p)  # np.diag(np.random.uniform(1 - spread_unif, 1 + spread_unif, p))
        # sigma0 = np.diag(sigma0)
        # sigma0[sigma0 == 0] = psi.numpy()
        # sigma0 = np.diag(sigma0)
        chol_params_0 = np.linalg.cholesky(sigma0).real
        is_neg_diag = np.where(np.diag(chol_params_0) < 0)
        if is_neg_diag[0].shape[0] >= 1:
            chol_params_0[is_neg_diag[0], :] *= -1
        np.fill_diagonal(chol_params_0, np.log(np.diag(chol_params_0)))
        idx_tril = np.tril_indices(p)
        chol_params_0 = chol_params_0[idx_tril]

        chol_params = t.tensor(chol_params_0).requires_grad_(True)
        mu_pep_mis = np.copy(mu0[idx_pep_mis])
        mu_pep_mis = t.tensor(mu_pep_mis).requires_grad_(True)
        n_pep_mis = mu_pep_mis.size(0)
        mu = t.tensor(mu0)

        if phi_known:
            phi = phi_init
            phi0 = phi0_init
            parameters = [mu_pep_mis, chol_params]
        else:
            log_phi = t.tensor([log_phi_init]).clone().requires_grad_(True)
            phi = t.exp(log_phi) + eps_phi
            phi0 = t.tensor([phi0_init]).clone().requires_grad_(True)
            parameters = [mu_pep_mis, chol_params, log_phi, phi0]

        optimizer = PyPirat.FullBatchLBFGS(parameters, line_search=line_search, lr=lr_init)
        optimizer.zero_grad()

        log_chol_m = get_matrix_from_tri_coefs(chol_params.detach().clone().numpy(),
                                               input_low=True, out_low=True, out_up=False)
        sigma_cur = get_sigma_from_log_chol(t.tensor(log_chol_m), eps_chol=eps_chol, eps_sig=eps_sig).numpy()
        if true_mu is not None:
            mu_errors = [np.sum((mu.detach().clone().numpy() - true_mu) ** 2)]
        if (true_sigma is not None) or (true_X is not None):
            if true_sigma is not None:
                sigma_errors = [np.sum((sigma_cur - true_sigma) ** 2)]
            if true_X is not None:
                ximp, var_imp = impute_from_params(X.numpy(), mu.detach().numpy(), sigma_cur,
                                                   phi.detach().clone().numpy(), phi0.detach().clone().numpy())
                imp_X_na = ximp[t.isnan(X).numpy() & ~np.isnan(true_X)]
                gt_X_na = true_X[t.isnan(X).numpy() & ~np.isnan(true_X)]
                ximp_errors = [np.sum((imp_X_na - gt_X_na) ** 2)]

        mu[idx_pep_mis] = mu_pep_mis

        obj = -compute_pen_llk(X, mu, chol_params, K, psi, phi0=phi0, phi=phi, eps_chol=eps_chol)
        objs = [-obj.item()]
        obj.backward()
        grad = [param.grad.flatten() for param in parameters]
        func_evals += 1

        x_old = [param.detach().clone().flatten() for param in parameters]
        f_old = obj
        # mu_list = [mu.detach().clone().numpy()]
        # sigma_list = [sigma_cur]
        for n_iter in range(maxiter):
            # print("####### ITER {} #######".format(n_iter))

            # define closure for line search
            if not phi_known:
                def closure():
                    optimizer.zero_grad()
                    phi = t.exp(log_phi) + eps_phi
                    mu = t.tensor(mu0)
                    mu[idx_pep_mis] = mu_pep_mis
                    loss_fn = -compute_pen_llk(X, mu, chol_params, K, psi, phi0=phi0, phi=phi, eps_chol=eps_chol,
                                               nsamples=nsamples)
                    return loss_fn
            else:
                def closure():
                    optimizer.zero_grad()
                    mu = t.tensor(mu0)
                    mu[idx_pep_mis] = mu_pep_mis
                    loss_fn = -compute_pen_llk(X, mu, chol_params, K, psi, phi0=phi0, phi=phi, eps_chol=eps_chol,
                                               nsamples=nsamples)
                    return loss_fn

            # perform line search step
            options = {'closure': closure, 'current_loss': obj, 'damping': damping, 'eta': eta, 'ls_debug': ls_debug,
                       'max_ls': max_ls, 'interpolate': interpolate, 'inplace': inplace}
            if line_search == 'Armijo':
                obj, lr, backtracks, clos_evals, desc_dir, fail = optimizer.step(options=options)

                # compute gradient at new iterate
                obj.backward()
                grad = optimizer._gather_flat_grad()

            elif line_search == 'Wolfe':
                obj, grad, lr, backtracks, clos_evals, grad_evals, desc_dir, fail = optimizer.step(options=options)

            mu = t.tensor(mu0)
            mu[idx_pep_mis] = mu_pep_mis

            log_chol_m = get_matrix_from_tri_coefs(chol_params.detach().clone().numpy(),
                                                   input_low=True, out_low=True, out_up=False)
            sigma_cur = get_sigma_from_log_chol(t.tensor(log_chol_m), eps_chol=eps_chol,
                                                eps_sig=eps_sig).numpy()
            if not phi_known:
                phi = t.exp(log_phi) + eps_phi # Update phi
            x_new = [param.detach().clone().flatten() for param in parameters]

            func_evals += clos_evals
            # compute quantities for checking convergence
            grads = [grad[:n_pep_mis], grad[n_pep_mis:]]
            grads_norm = t.max(t.tensor([t.norm(gradi) for gradi in grads]))
            x_dist = t.max(t.tensor([
                t.norm(x_new[i] - x_old[i])/t.norm(x_old[i]) for i in range(len(x_new))
            ]))
            # xdiffs.append(x_dist.item())

            f_dist = t.abs(obj - f_old) / t.abs(f_old)

            # print data
            if out:
                print('  %.3e  |  %.3e  |  %.3e  |  %.3e  |  %.3e  |  %.3e  ' % (
                n_iter + 1, obj.item(), grads_norm.item(), x_dist.item(), clos_evals, lr))

            objs.append(-obj.item())

            if true_mu is not None:
                mu_errors.append(np.sum((mu.detach().clone().numpy() - true_mu) ** 2))
            if (true_sigma is not None) or (true_X is not None):
                if true_sigma is not None:
                    sigma_errors.append(np.sum((sigma_cur - true_sigma) ** 2))
                if true_X is not None and n_iter % 10 == 0:
                    ximp, var_imp = impute_from_params(
                        X.numpy(), mu.detach().numpy(), sigma_cur, phi.detach().clone().numpy(), phi0.detach().clone().numpy())
                    imp_X_na = ximp[t.isnan(X).numpy() & ~np.isnan(true_X)]
                    gt_X_na = true_X[t.isnan(X).numpy() & ~np.isnan(true_X)]
                    ximp_errors.append(np.sum((imp_X_na - gt_X_na) ** 2))

            # mu_list.append(mu.detach().clone().numpy())
            # sigma_list.append(sigma_cur)

            # stopping criterion
            if (t.isnan(obj) or n_iter == maxiter - 1) or obj <= -1e8 or fail:
                TRUEFAIL = True
                if obj <= -1e8:
                    error_msg = "ill_posed"
                elif t.isnan(obj):
                    error_msg = "loss_eq_nan"
                elif fail:
                    error_msg = "ls_failed"
                else:
                    error_msg = "max_iter"
                reason = None
                break
            elif grads_norm.item() < tol_grad or x_dist < tol_param or f_dist < tol_obj:
                if f_dist < tol_obj:
                    reason = "obj"
                if x_dist < tol_param:
                    reason = "param"
                if grads_norm.item() < tol_grad:
                    reason = "grad"
                if np.nanmax(sigma_cur) >= 200:
                    error_msg = "var_exploded"
                else:
                    error_msg = "success"
                break

            x_old = [param.detach().clone().flatten() for param in parameters]
            f_old.copy_(obj)
        lr_init = lr_init/2.

    # print("NB TRY :", nb_try)

    log_chol_m = get_matrix_from_tri_coefs(chol_params.detach().clone().numpy(),
                                           input_low=True, out_low=True, out_up=False)
    sigma_cur = get_sigma_from_log_chol(t.tensor(log_chol_m), eps_chol=eps_chol, eps_sig=eps_sig)
    ximp, var_imp = impute_from_params(X.numpy(), mu.detach().numpy(), sigma_cur.detach().numpy(), phi.detach().clone().numpy(),
                                       phi0.detach().clone().numpy())
    if true_X is not None:
        imp_X_na = ximp[t.isnan(X).numpy() & ~np.isnan(true_X)]
        gt_X_na = true_X[t.isnan(X).numpy() & ~np.isnan(true_X)]
        ximp_errors.append(np.sum((imp_X_na - gt_X_na) ** 2))


    time_cur = time.process_time() - time_cur

    # print summary
    if verbose:
        print('==================================== Summary ======================================')
        if TRUEFAIL:
            print("------- FAIL ----------")
            print(error_msg)
        else:
            print("------- SUCCESS ----------")
        print('Iterations:', n_iter + 1)
        print('Function Evaluations:', func_evals)
        print('Time:', time_cur)
        print('F:', obj.item())
        print('x_dist', x_dist)
        print('||g||:', t.norm(grad).item())
        print('===================================================================================')

    # res_dic = {'Xhat': ximp, 'varXhat': var_imp, 'phi0': phi0.detach().numpy(), 'phi': phi.detach().numpy(),
    #            'llks': np.array(objs), 'mu': np.squeeze(np.array(mu_list))[:,idx_pep_mis],
    #            'sigma': np.squeeze(np.array(sigma_list)),
    #            'error_msg': error_msg, "sd_obj": np.sqrt(np.var(lobj)), 'reason': reason, 'xdiffs': xdiffs,
    #            'obj_comp': np.array(lobj)}

    res_dic = {'Xhat': ximp, 'varXhat': var_imp, 'phi0': phi0.detach().numpy(), 'phi': phi.detach().numpy(),
               'llks': np.array(objs), 'mu': mu.detach().numpy(),
               'sigma': sigma_cur.detach().numpy(),
               'error_msg': error_msg, 'reason': reason} # 'objs_end': np.array([lobj])}

    if true_mu is not None:
        res_dic['mu_errors'] = np.array(mu_errors)
    if true_sigma is not None:
        res_dic['sigma_errors'] = np.array(sigma_errors)
    if true_X is not None:
        res_dic['ximp_errors'] = np.array(ximp_errors)
        res_dic['Xgt'] = true_X

    return res_dic
