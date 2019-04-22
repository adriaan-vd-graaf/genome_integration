import scipy.stats
import numpy as np
import statsmodels.api  as sm
from sklearn.linear_model import BayesianRidge, LassoCV


def remove_highly_correlated_snps(r_sq_mat, r_sq_threshold=0.95):
    """

    This iteratively selects variants that are less correlated to each other than r_sq_threshold
    removes SNPs that are more correlated that `r_sq_threshold`.
    First SNP in r_sq_mat is always retained.


    :param r_sq_mat: squared pearson correlation matrix of variants, shape (m x m)
    :param r_sq_threshold: threshold for `r_sq_mat`, float in [0,1]
    :return:  boolean array of shape (m) indicating which SNPs are lower in LD than `r_sq_threshold`
    """

    tril = np.tril(r_sq_mat, -1)
    mask = np.ones((r_sq_mat.shape[0]), dtype=bool)
    # I believe this is allowed because we are also using permutations sometimes.
    i = 0
    while i < r_sq_mat.shape[0]:
        if mask[i]:
            mask[tril[:,i] > r_sq_threshold] = False
            i += 1
        else:
            i += 1

    return mask


def mask_instruments_in_ld(r_sq_mat,
                           instruments,
                           upper_r_sq_thresh=0.99,
                           lower_r_sq_thresh=0.1,
                           prune_r_sq_thresh=0.95,
                           shuffle_positions=False):
    """


    :param r_sq_mat: squared pearson correlation matrix of variants, shape (m x m)
    :param instruments: instruments (variants) used in MR
    :param upper_r_sq_thresh: maximum correlation from instrument threshold for `r_sq_mat`, float in [0,1]
    :param lower_r_sq_thresh: minimum from instrument threshold for `r_sq_mat`, float in [0,1]
    :param prune_r_sq_thresh: threshold from which to remove highly correlated SNPs.
    :return:
    """


    masked_instruments = np.zeros((r_sq_mat.shape[0]), dtype=bool)

    #if there are any NA's in the r_sq matrix
    masked_instruments[np.any(np.isnan(r_sq_mat), axis=0)] = True

    for instrument in instruments:
        tmp_mask = r_sq_mat[instrument,:] < upper_r_sq_thresh
        tmp_mask = np.logical_and(tmp_mask, r_sq_mat[instrument,:] > lower_r_sq_thresh)

        masked_instruments[tmp_mask] = True


    if shuffle_positions:
        int_indices = np.where(masked_instruments)[0]
        np.random.shuffle(int_indices)
        highly_correlated = remove_highly_correlated_snps(r_sq_mat[int_indices, :][:, int_indices],
                                                          prune_r_sq_thresh)

    else:
        highly_correlated = remove_highly_correlated_snps(r_sq_mat[masked_instruments,:][:,masked_instruments],
                                                          prune_r_sq_thresh)

    masked_instruments[masked_instruments] = highly_correlated

    return masked_instruments


def mr_link_ols(outcome_geno,
                r_sq_mat,
                exposure_betas,
                causal_exposure_indices,
                outcome_phenotype,
                upper_r_sq_threshold=0.99
                ):

    masked_instruments = mask_instruments_in_ld(r_sq_mat,
                                                causal_exposure_indices,
                                                upper_r_sq_threshold)

    geno_masked = outcome_geno[masked_instruments,:].transpose()

    design_mat = np.zeros( (geno_masked.shape[0], geno_masked.shape[1]+1), dtype=float)
    design_mat[:,0] = (outcome_geno[causal_exposure_indices,:].transpose() @  exposure_betas) / exposure_betas.shape[0]
    design_mat[:,np.arange(1, design_mat.shape[1])] = geno_masked / np.sqrt(geno_masked.shape[1])

    ols_fit = sm.OLS(endog=outcome_phenotype, exog=design_mat).fit()

    return ols_fit.params[0], ols_fit.bse[0], ols_fit.pvalues[0]


def mr_link_ridge(outcome_geno,
                 r_sq_mat,
                 exposure_betas,
                 causal_exposure_indices,
                 outcome_phenotype,
                 upper_r_sq_threshold=0.99,
                 use_hat_matrix_for_p_determination=False
                 ):

    masked_instruments = mask_instruments_in_ld(r_sq_mat,
                                                causal_exposure_indices,
                                                upper_r_sq_threshold)


    geno_masked = outcome_geno[masked_instruments,:].transpose()

    ridge_fit = BayesianRidge(fit_intercept=False)

    design_mat = np.zeros( (geno_masked.shape[0], geno_masked.shape[1]+1), dtype=float)
    design_mat[:,0] = (outcome_geno[causal_exposure_indices,:].transpose() @ exposure_betas) / exposure_betas.shape[0]
    design_mat[:,np.arange(1, design_mat.shape[1])] = geno_masked / np.sqrt(geno_masked.shape[1])

    ridge_fit.fit(design_mat, outcome_phenotype)

    t_stat = np.abs(ridge_fit.coef_[0] / np.sqrt(ridge_fit.sigma_[0, 0]))

    if not use_hat_matrix_for_p_determination:

        p_val = 2 * scipy.stats.norm.sf(t_stat)

    else:
        # this produces a very big speed decrease, but should be done when the design matrix
        # contains more columns than there are individuals in the cohort.
        hat_mat = design_mat @ np.linalg.inv(ridge_fit.sigma_) @ design_mat.T
        deg_freedom = design_mat.shape[0] - np.trace(2*hat_mat - hat_mat @ hat_mat.T)
        print(deg_freedom)
        p_val = 2* scipy.stats.t.sf(t_stat, deg_freedom)

    return ridge_fit.coef_[0], np.sqrt(ridge_fit.sigma_[0,0]), p_val


def mr_link_ridge_resample_tags(outcome_geno,
                                r_sq_mat,
                                exposure_betas,
                                causal_exposure_indices,
                                outcome_phenotype,
                                upper_r_sq_threshold=0.99,
                                n_resamples = 100
                                ):

    masked_instruments = mask_instruments_in_ld(r_sq_mat,
                                                causal_exposure_indices,
                                                upper_r_sq_threshold)


    ridge_fit = BayesianRidge(fit_intercept=False)
    geno_masked = outcome_geno[masked_instruments,:].transpose()


    design_mat = np.zeros( (geno_masked.shape[0], geno_masked.shape[1]+1), dtype=float)
    design_mat[:,0] = (outcome_geno[causal_exposure_indices,:].transpose() @ exposure_betas) / exposure_betas.shape[0]
    design_mat[:,np.arange(1, design_mat.shape[1])] = geno_masked / np.sqrt(geno_masked.shape[1])

    ridge_fit.fit(design_mat, outcome_phenotype)

    original_z_stat = ridge_fit.coef_[0] / np.sqrt(ridge_fit.sigma_[0, 0])
    p_val = 2 * scipy.stats.norm.sf(np.abs(original_z_stat))

    original_regression = ridge_fit.coef_[0], np.sqrt(ridge_fit.sigma_[0,0]), p_val

    all_zs = np.zeros(n_resamples, dtype=float)
    i = 0
    faults = 0
    while i < n_resamples:
        try:
            sampled_tags = mask_instruments_in_ld(r_sq_mat,
                                                  causal_exposure_indices,
                                                  upper_r_sq_threshold,
                                                  shuffle_positions=True)

            geno_masked = outcome_geno[sampled_tags, :].transpose()
            design_mat = np.zeros((geno_masked.shape[0], geno_masked.shape[1] + 1), dtype=float)
            design_mat[:, 0] = (outcome_geno[causal_exposure_indices, :].transpose() @ exposure_betas) / \
                               exposure_betas.shape[0]
            design_mat[:, np.arange(1, design_mat.shape[1])] = geno_masked / np.sqrt(geno_masked.shape[1])

            fit = ridge_fit.fit(design_mat, outcome_phenotype)
            sampled_z_stat = ridge_fit.coef_[0] / np.sqrt(ridge_fit.sigma_[0, 0])
            all_zs[i] = sampled_z_stat
            i += 1

        except np.linalg.LinAlgError:
            faults += 1
            if faults >= (0.1*n_resamples):
                raise ValueError(f"Did ridge tag SNP resampling, but SVD did not "
                                 f"converge in more than {faults} out of {i} estimations")


    new_se = np.std(all_zs)
    new_mean = np.mean(all_zs)
    direction_concordance = np.sum(np.sign(all_zs) != np.sign(original_z_stat)) / n_resamples

    new_p = 2 * scipy.stats.t.sf(np.abs(new_mean), n_resamples-1)
    return ((new_mean, new_se, new_p),
            (original_regression[0], original_regression[1], direction_concordance))


def mr_link_ridge_cv(outcome_geno,
                     r_sq_mat,
                     exposure_betas,
                     causal_exposure_indices,
                     outcome_phenotype,
                     upper_r_sq_threshold=0.99
                     ):

    masked_instruments = mask_instruments_in_ld(r_sq_mat,
                                                causal_exposure_indices,
                                                upper_r_sq_threshold)


    geno_masked = outcome_geno[masked_instruments,:].transpose()

    bayesian_ridge_fit = BayesianRidge(fit_intercept=False)

    design_mat = np.zeros( (geno_masked.shape[0], geno_masked.shape[1]+1), dtype=float)
    design_mat[:,0] = (outcome_geno[causal_exposure_indices,:].transpose() @ exposure_betas) / exposure_betas.shape[0]
    design_mat[:,np.arange(1, design_mat.shape[1])] = geno_masked / np.sqrt(geno_masked.shape[1])

    full_fit = bayesian_ridge_fit.fit(X=design_mat, y=outcome_phenotype)

    ridge_fit = BayesianRidge(fit_intercept=False)

    n_splits = 100
    save_array = np.zeros(shape=(n_splits, 3), dtype=float)

    for i in range(n_splits):
        train_index = np.random.choice(design_mat.shape[0],int(np.floor(design_mat.shape[0] * .5)), replace=True)
        fit = ridge_fit.fit(X=design_mat[train_index,:], y=outcome_phenotype[train_index])
        save_array[i,:] = fit.coef_[0] / np.sqrt(fit.sigma_[0,0]), fit.coef_[0], None


    mean_t = full_fit.coef_[0] / np.sqrt(full_fit.sigma_[0,0])
    se_t = np.std(save_array[:,0], ddof=1)
    if np.mean(save_array[:, 0]) > 0.0:

        quantile_p_val_t = (np.sum(save_array[:,0] < 0) / n_splits) * 2
    else:
        quantile_p_val_t = (np.sum(save_array[:, 0] > 0) / n_splits) * 2

    mean_b = full_fit.coef_[0]
    se_b = np.std(save_array[:, 1], ddof=1)
    if np.mean(save_array[:, 1]) > 0:
        quantile_p_val_b = (np.sum(save_array[:, 1] < 0) / n_splits) * 2
    else:
        quantile_p_val_b = (np.sum(save_array[:, 1] > 0) / n_splits) * 2

    return [[mean_t, se_t, quantile_p_val_t], [mean_b, se_b,quantile_p_val_b]]


def mr_link_bootstrapped(outcome_geno,
                         r_sq_mat,
                         exposure_betas,
                         causal_exposure_indices,
                         outcome_phenotype,
                         upper_r_sq_threshold=0.99,
                         bootstraps = 50
                         ):


    masked_instruments = mask_instruments_in_ld(r_sq_mat,
                                                causal_exposure_indices,
                                                upper_r_sq_threshold)


    geno_masked = outcome_geno[masked_instruments,:].transpose()

    ridge_fit = BayesianRidge(fit_intercept=False)

    design_mat = np.zeros( (geno_masked.shape[0], geno_masked.shape[1]+1), dtype=float)
    design_mat[:,0] = (outcome_geno[causal_exposure_indices,:].transpose() @ exposure_betas) / exposure_betas.shape[0]
    design_mat[:,np.arange(1, design_mat.shape[1])] = geno_masked / np.sqrt(geno_masked.shape[1])

    ridge_fit.fit(design_mat, outcome_phenotype)

    t_stat = np.abs(ridge_fit.coef_[0] / np.sqrt(ridge_fit.sigma_[0, 0]))

    p_val = 2 * scipy.stats.norm.sf(t_stat)

    if type(bootstraps) is int:
        max_iter = bootstraps
        checkpoints = {bootstraps}
    else:
        max_iter = max(bootstraps)
        checkpoints = set(bootstraps)

    bootstrapping_results = []
    returns = {}

    for i in range(max_iter):
        selected_individuals = np.random.choice(outcome_phenotype.shape[0],
                                                int(outcome_phenotype.shape[0]),
                                                replace=True
                                                )

        bootstrapped_fit = ridge_fit.fit(design_mat[selected_individuals,:], outcome_phenotype[selected_individuals])
        bootstrapping_results += [bootstrapped_fit.coef_[0]]

        if i+1 in checkpoints:
            obs_se = np.std(bootstrapping_results, ddof=1)
            p_val = 2 * scipy.stats.norm.sf(np.abs(ridge_fit.coef_[0]) / obs_se)

            returns[i+1] = (ridge_fit.coef_[0], obs_se, p_val)

    print(np.quantile(bootstrapping_results, [0.025, 0.975]))

    if type(bootstraps) is int:
        return returns[bootstraps]
    else:
        return returns


def mr_link_lasso(outcome_geno,
                  r_sq_mat,
                  exposure_betas,
                  causal_exposure_indices,
                  outcome_phenotype,
                  upper_r_sq_threshold = 0.95,
                  ):

    masked_instruments = mask_instruments_in_ld(r_sq_mat,
                                                causal_exposure_indices,
                                                upper_r_sq_threshold,
                                                lower_r_sq_thresh=0.1,
                                                prune_r_sq_thresh=0.95)

    lasso_fit = LassoCV(fit_intercept=False, cv=5, max_iter=5000, eps=0.001)
    masked_geno = outcome_geno[masked_instruments,:]
    lasso_fit.fit(X=masked_geno.T, y=outcome_phenotype)

    to_keep = lasso_fit.coef_ != 0.0

    print("check3")
    new_mat = np.zeros((masked_geno.shape[1], np.sum(to_keep)+1), dtype=float)
    new_mat[:,0] = ((outcome_geno[causal_exposure_indices, :].transpose() @ exposure_betas) / exposure_betas.shape[0])
    new_mat[:,np.arange(1, np.sum(to_keep)+1)] = masked_geno[to_keep,:].transpose()

    print(sum(to_keep))

    print("check4")

    ols_fit = sm.OLS(endog=outcome_phenotype, exog=new_mat).fit()

    return ols_fit.params[0], ols_fit.bse[0], ols_fit.pvalues[0]
