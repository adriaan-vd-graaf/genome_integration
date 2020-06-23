import copy
import scipy.stats
import numpy as np
import statsmodels.api  as sm
from sklearn.linear_model import BayesianRidge

"""
 Unused functionality was deleted on the 17th of June, if you are interested in it, 
 please go back to a commit before that date.
"""

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
                                  upper_r_sq_threshold=0.99,
                                  output_selected_variants = False):
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

    if output_selected_variants:
        return design_mat, masked_instruments
    else:
        return design_mat


def mr_link_ols(outcome_geno,
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
    Please note that the p value and se is uncorrected. so these are usually _very_ conservative.
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


def make_knockoff_genotypes(plinkfile, tmp_loc='tmp_files_fastphase'):
    """
    This makes knockoff genotypes.
    in the first step, fastphase is run
    in the second step, the knockoff methodology is performed.

    it will return a genotype matrix e

    :param plinkfile: PlinkFile object from which knockoffs will be created.
    :param tmp_loc:  temporary location where to store run files/
    :return: a genotype matrix of the same size as the PlinkFile.maf object
    """

    fastphase_file = f'{tmp_loc}_input.inp'
    with open(fastphase_file, 'w') as f:
        f.write(f'{len(plinkfile.fam_data.sample_names)}\n')
        f.write(f'{len(plinkfile.bim_data.snp_names)}\n')
        position_string = \
            f'P {" ".join(plinkfile.bim_data.bim_results[x].position for x in plinkfile.bim_data.snp_names)}'
        f.write(f'{position_string}\n')
        for i, sample in enumerate(plinkfile.fam_data.sample_names):
            f.write(f'# {sample}\n')
            genotypes = plinkfile.genotypes[i, :]
            first_line = ''.join(['0' if int(x) == 0 else '1' for x in genotypes])
            second_line = ''.join(['0' if int(x) == 2 else '1' for x in genotypes])
            f.write(f'{first_line}\n')
            f.write(f'{second_line}\n')

    fastphase_out = f'{tmp_loc}_fastphase_output'

    subprocess.run(['fastPHASE',
                    '-Pp',
                    '-T1',
                    '-K15',  # this is the amount of haplotype clusters there are.
                    '-g',
                    '-H-4',
                    '-C20',  # this is the number of iterations
                    f'-o{fastphase_out}',
                    fastphase_file
                    ], check=True)

    r_file = f"{fastphase_out}_rhat.txt"
    alpha_file = f"{fastphase_out}_alphahat.txt"
    theta_file = f"{fastphase_out}_thetahat.txt"
    char_file = f"{fastphase_out}_origchars"
    hmm = loadHMM(r_file, alpha_file, theta_file, char_file, compact=True)

    knockoffs_gen = SNPknock.knockoffGenotypes(hmm['r'], hmm['alpha'], hmm['theta'])
    Xk = knockoffs_gen.sample(plinkfile.genotypes)

    # clean up.
    [os.remove(x) for x in [r_file, alpha_file, theta_file, char_file, fastphase_file]]
    return Xk


def knockoff_filter_threshold(w_vector, fdr=0.05, offset=1):
    """
    R function this is based off. from the SNPknock package by Matteo Sesia

    knockoff.threshold < - function(W, fdr=0.10, offset=1)
    {
    if (offset != 1 & & offset != 0)
    {
        stop('Input offset must be either 0 or 1')
    }
    ts = sort(c(0, abs(W)))
    ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
    ok = which(ratio <= fdr)
    ifelse(length(ok) > 0, ts[ok[1]], Inf)
    }

    :param w_vector:
    :param fdr:
    :param offset:
    :return: threshold for the W vector.
    """

    if offset not in [0, 1]:
        raise ValueError("Offset should be in the set 0 or 1")
    ts = np.asarray(np.abs([0] + list(w_vector)), dtype=float)
    ratios = np.asarray([(offset + np.sum(w_vector <= -t)) / max([1, np.sum(w_vector >= t)]) for t in ts], dtype=float)
    threshold = np.min(ts[np.logical_and(ratios < fdr, ts > 0)])
    return threshold


def mr_link_knockoffs(
        outcome_plinkfile,
        scaled_outcome_geno,
        outcome_ld_r_sq,
        beta_effects,
        iv_selection,
        outcome_phenotypes,
        tmp_file='tmp_for_mr_link_knockoffs',
        upper_r_sq_threshold=0.99
):
    design_matrix, tag_indices = make_mr_link_design_matrix(scaled_outcome_geno,
                                                             outcome_ld_r_sq,
                                                             beta_effects,
                                                             iv_selection,
                                                             upper_r_sq_threshold,
                                                             output_selected_variants=True)

    # all indices contains the tag snps and the ivs.
    tag_snp_names = np.asarray(outcome_plinkfile.bim_data.snp_names, dtype=str)[tag_indices]
    iv_snp_names = np.asarray(outcome_plinkfile.bim_data.snp_names, dtype=str)[iv_selection]

    snps_to_keep = list(iv_snp_names) + list(tag_snp_names)

    pruned_outcome_plinkfile = copy.copy(outcome_plinkfile)
    pruned_outcome_plinkfile.prune_for_a_list_of_snps(snps_to_keep)
    knockoff_genotypes = make_knockoff_genotypes(pruned_outcome_plinkfile, tmp_loc=tmp_file + "_knockoff_intermediates")

    iv_indices = np.asarray([pruned_outcome_plinkfile.bim_data.snp_names.index(x) for x in iv_snp_names], dtype=int)
    tag_indices = np.asarray([pruned_outcome_plinkfile.bim_data.snp_names.index(x) for x in tag_snp_names], dtype=int)

    knockoff_design_mat = np.zeros(design_matrix.shape, dtype=float)
    knockoff_design_mat[:, 0] = (knockoff_genotypes[:, iv_indices] @ beta_effects) / beta_effects.shape[0]
    knockoff_design_mat[:, np.arange(1, knockoff_design_mat.shape[1])] = knockoff_genotypes[:, tag_indices] / np.sqrt(
        tag_indices.shape[0])

    mr_link_orig_and_knockoff_joint = BayesianRidge(fit_intercept=False)
    mr_link_orig_and_knockoff_joint.fit(np.concatenate([design_matrix, knockoff_design_mat], axis=1),
                                        outcome_phenotypes)

    all_t_stats = np.abs(
        mr_link_orig_and_knockoff_joint.coef_ / np.sqrt(np.diag(mr_link_orig_and_knockoff_joint.sigma_)))
    # all_ws = all_t_stats[:design_matrix.shape[1]] - all_t_stats[design_matrix.shape[1]:] #W in the Candes / Barbes nomenclature
    all_p_vals = 2 * scipy.stats.norm.sf(all_t_stats)

    print(f'mr_link orig joint {mr_link_orig_and_knockoff_joint.coef_[0]:.3f}, {all_p_vals[0]:.3e}')
    print(
        f'mr_link knockoffs joint {mr_link_orig_and_knockoff_joint.coef_[design_matrix.shape[1]]:.3f}, {all_p_vals[design_matrix.shape[1]]:.3e}')

    all_ws = np.abs(mr_link_orig_and_knockoff_joint.coef_[:144]) - np.abs(mr_link_orig_and_knockoff_joint.coef_[144:])
    w_threshold = knockoff_filter_threshold(all_ws, fdr=0.05, offset=1)

    return mr_link_orig_and_knockoff_joint, all_ws, w_threshold, pruned_outcome_plinkfile.bim_data.snp_names