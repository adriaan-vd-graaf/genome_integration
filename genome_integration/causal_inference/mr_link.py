
import scipy.stats
import numpy as np
import statsmodels.api  as sm
from sklearn.linear_model import BayesianRidge


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
    i = 0
    while i < r_sq_mat.shape[0]:
        if mask[i]:
            mask[tril[:,i] > r_sq_threshold] = False
            i += 1
        else:
            i += 1

    return mask


def mask_instruments_in_ld(r_sq_mat, instruments, upper_r_sq_thresh=0.99, lower_r_sq_thresh=0.1, prune_r_sq_thresh=0.95):
    """


    :param r_sq_mat: squared pearson correlation matrix of variants, shape (m x m)
    :param instruments: instruments (variants) used in MR
    :param upper_r_sq_thresh: maximum correlation from instrument threshold for `r_sq_mat`, float in [0,1]
    :param lower_r_sq_thresh: minimum from instrument threshold for `r_sq_mat`, float in [0,1]
    :param prune_r_sq_thresh: threshold from which to remove highly correlated SNPs.
    :return:
    """


    masked_instruments = np.zeros((r_sq_mat.shape[0]), dtype=bool)
    for instrument in instruments:
        tmp_mask = r_sq_mat[instrument,:] < upper_r_sq_thresh
        tmp_mask = np.logical_and(tmp_mask, r_sq_mat[instrument,:] > lower_r_sq_thresh)

        masked_instruments[tmp_mask] = True

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
    p_val = 2 * scipy.stats.norm.sf(np.abs(ridge_fit.coef_[0] / np.sqrt(ridge_fit.sigma_[0, 0])))

    return ridge_fit.coef_[0], np.sqrt(ridge_fit.sigma_[0,0]), p_val
