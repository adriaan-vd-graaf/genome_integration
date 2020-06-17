import scipy.stats
import numpy as np
import statsmodels.api  as sm
from sklearn.linear_model import BayesianRidge, LassoCV, LogisticRegression


def remove_highly_correlated_snps(r_sq_mat, r_sq_threshold=0.95):
    """
    This iteratively selects variants that are less correlated to each other than r_sq_threshold
    removes SNPs that are more correlated that `r_sq_threshold`.

    First SNP in r_sq_mat is always retained, later ones are pruned.
    Of note, in simulations and in real data, when permuting the ordering of r_sq_mat, we found no meaningful
    differences in MR-link results. So

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


def mask_instruments_in_ld(r_sq_mat,
                           instruments,
                           upper_r_sq_thresh=0.99,
                           lower_r_sq_thresh=0.1,
                           prune_r_sq_thresh=0.95,
                           shuffle_positions=False):
    """
    Masks instruments that are in LD [upper_r_sq_threshold, lower_r_sq_threshold] with instrumental variables of the
    exposure. As well as being in high LD (> prune_r_sq_threshold) with itself.

    :param r_sq_mat: squared pearson correlation matrix of variants, shape (m x m)
    :param instruments: indices () of the instruments (variants) used for the exposure.
    :param upper_r_sq_thresh: maximum correlation from instrument threshold for `r_sq_mat`, float in [0,1]
    :param lower_r_sq_thresh: minimum from instrument threshold for `r_sq_mat`, float in [0,1]
    :param prune_r_sq_thresh: threshold from which to remove highly correlated SNPs. float in [0,1]

    :return: Returns a logical vector of length m representing variants that are in high to low LD with the
    instrumental variables and in high LD with other variants.
    """

    masked_variants = np.zeros((r_sq_mat.shape[0]), dtype=bool)

    #if there are any NA's in the r_sq matrix
    masked_variants[np.any(np.isnan(r_sq_mat), axis=0)] = True

    for instrument in instruments:
        tmp_mask = r_sq_mat[instrument,:] < upper_r_sq_thresh
        tmp_mask = np.logical_and(tmp_mask, r_sq_mat[instrument,:] > lower_r_sq_thresh)

        masked_variants[tmp_mask] = True

    # This shuffles the LD matrix, just to be sure.
    # Analysis of shuffling the ordering of pruning, did not to change the results in preliminary analysis
    # Therefore it is not used.
    if shuffle_positions:
        int_indices = np.where(masked_variants)[0]
        np.random.shuffle(int_indices)
        highly_correlated = remove_highly_correlated_snps(r_sq_mat[int_indices, :][:, int_indices],
                                                          prune_r_sq_thresh)

    else:
        highly_correlated = remove_highly_correlated_snps(r_sq_mat[masked_variants,:][:,masked_variants],
                                                          prune_r_sq_thresh)

    masked_variants[masked_variants] = highly_correlated

    return masked_variants


def make_mr_link_design_matrix(outcome_geno,
                                  r_sq_mat,
                                  exposure_betas,
                                  causal_exposure_indices,
                                  upper_r_sq_threshold=0.99):
    """

    :param outcome_geno: genotype matrix of the outcome
    :param r_sq_mat: R^2 matrix of all the genotypes of the outcome
    :param exposure_betas: beta estimates of the exposure instrumental variables
    :param causal_exposure_indices: indices of the exposure instrumental variables
    :param upper_r_sq_threshold: the upper r_sq threshold for which the variants around the IVs are pruned
    :return: a design matrix for use in MR-link.
    """

    masked_instruments = mask_instruments_in_ld(r_sq_mat,
                                                causal_exposure_indices,
                                                upper_r_sq_threshold)


    geno_masked = outcome_geno[:,masked_instruments]

    design_mat = np.zeros( (geno_masked.shape[0], geno_masked.shape[1]+1), dtype=float)
    design_mat[:,0] = (outcome_geno[:,causal_exposure_indices] @ exposure_betas) / exposure_betas.shape[0]
    design_mat[:,np.arange(1, design_mat.shape[1])] = geno_masked / np.sqrt(geno_masked.shape[1])

    return design_mat


def mr_link_ols_old(outcome_geno,
                r_sq_mat,
                exposure_betas,
                causal_exposure_indices,
                outcome_phenotype,
                upper_r_sq_threshold=0.99
                ):
    """
    Does MR-link solved by ordinary least squares.

    :param outcome_geno: outcome genotypes
    :param r_sq_mat: R^2 matrix in order of genotypes of outcome geno
    :param exposure_betas: beta estimates of the exposure instrumental variables.
    :param causal_exposure_indices:  indices of the exposure instrumental variables.
    :param outcome_phenotype: outcome phenotype vector
    :param upper_r_sq_threshold: the upper r_sq threshold for which the variants around the IVs are pruned.
    :return: beta, se and p value estimate of the MR-link estimate
    """

    design_mat = make_mr_link_design_matrix(outcome_geno,
                                            r_sq_mat,
                                            exposure_betas,
                                            causal_exposure_indices,
                                            upper_r_sq_threshold)

    ols_fit = sm.OLS(endog=outcome_phenotype, exog=design_mat).fit()

    return ols_fit.params[0], ols_fit.bse[0], ols_fit.pvalues[0]


def mr_link_ridge(outcome_geno,
                 r_sq_mat,
                 exposure_betas,
                 causal_exposure_indices,
                 outcome_phenotype,
                 upper_r_sq_threshold=0.99
                 ):
    """

    Does MR-link solved by ridge regression.
    Please note that the p value and se is uncorrected. so these is usually _very_ conservative.
    See the MR-link manuscript for details.

    :param outcome_geno: outcome genotypes
    :param r_sq_mat: R^2 matrix in order of genotypes of outcome geno
    :param exposure_betas: beta estimates of the exposure instrumental variables.
    :param causal_exposure_indices:  indices of the exposure instrumental variables.
    :param outcome_phenotype: outcome phenotype vector
    :param upper_r_sq_threshold: the upper r_sq threshold for which the variants around the IVs are pruned.
    :return: beta, se and p value estimate of the MR-link estimate
    """

    design_mat = make_mr_link_design_matrix(outcome_geno,
                                            r_sq_mat,
                                            exposure_betas,
                                            causal_exposure_indices,
                                            upper_r_sq_threshold)

    ridge_fit = BayesianRidge(fit_intercept=False)

    ridge_fit.fit(design_mat, outcome_phenotype)

    t_stat = np.abs(ridge_fit.coef_[0] / np.sqrt(ridge_fit.sigma_[0, 0]))

    p_val = 2 * scipy.stats.norm.sf(t_stat)

    return ridge_fit.coef_[0], np.sqrt(ridge_fit.sigma_[0,0]), p_val