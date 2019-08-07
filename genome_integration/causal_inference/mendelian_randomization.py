"""
This class is used to do a 'simple'
inverse variance estimate on independent effects.
"""

__author__      = "Adriaan van der Graaf"

import math
import copy
import warnings
import numpy as np
import scipy.stats
from statsmodels.regression.linear_model import WLS
from statsmodels.tools.tools import add_constant


from .. import utils
from .. import association


class MendelianRandomization:
    """
    This class is a base class of most summary statistic based MR analyses.


    Attributes
    ----------


    estimation_done: bool
        boolean indicating if an estimation was done.

    estimation_data: empty list
        initialized as an empty list,
        but as estimates are added, will contain tuples of length 2 with a beta and a standard error.

    estimation_snps = : empty list
        initialized as an empty list, but will be filled with snp information as estimates are added. can be left empty.

    outcome_tuples: empty list
        initialized as an empty list, can be filled with the summary statistics of the outcome
        in the same way as estimation data

    exposure_tuples: empty list
        initialized as an empty list, can be filled with the summary statistics of the exposure
        in the same way as estimation data

    q_test_indices_remaining: list of ints
        initialized as empty, but as the q test removes estimates, it shows which estimates are removed.

    _ivw_intermediate_top: list of floats
        this is an internal variable that should not be used.

    _ivw_intermediate_bottom list of floats
        this is an internal variable that should not be used.


    Methods
    -------
    add_estimate(self, beta_se_tuple, variant_name, pos, chr)
        adds a single estimate for MR estimation.
        beta_se_tuple is a tuple of length 2, containing the beta estimate and se(beta).
        The variant name, pos and chr are information that denote where the variant comes from.

    do_ivw_estimation(self)
        Does the IVW estimation based on the single estimates in estimation data.

    do_ivw_estimation_on_estimate_vector(self, estimation_vec)
        Does an IVW estimate on an estimate vector, without storing any information in the class.
        estimation_vec is a list of tuples with beta and se per variant that's being used.

    get_ivw_estimates(self)
        Deprecated, mirrors the do_ivw_estimation() method.

    do_smr_estimate(self, exposure_tuple, outcome_tuple)
        Does an SMR estimate of causality. SMR uses the variance of the outcome summary statistic and the variance of
        the exposure summary statistic to estimate significance.
        exposure_tuple is a beta se tuple of the exposure summary statistics
        outcome_tuple is a beta se tuple of the outcome summary statistics

    do_smr_estimate(self, exposure_tuple, outcome_tuple)
        Does an SMR estimate of causality. SMR uses the variance of the outcome summary statistic and the variance of
        the exposure summary statistic to estimate significance.
        exposure_tuple is a beta se tuple of the exposure summary statistics
        outcome_tuple is a beta se tuple of the outcome summary statistics

    do_and_add_smr_estimation(self, exposure_tuple, outcome_tuple, variant_name, pos, chr)
        Adds an SMR estimate from the exposure and outcome tuples (beta, se).

    do_chochrans_q_meta_analysis(self, p_value_threshold):
        Does a chochrans q meta analysis of the results at a certain p value threhold.

    do_egger_regression(self):
        Does the original Egger regression with a single variance term. (weighted by the outcome se)

    do_egger_regression_single_variance_term(self):
        Does egger regression with a single variance term. (weighted only by the outcome se)

    do_egger_regression_two_variance_term(self)
        Does egger regression with a double variance term. (weighted by both the outcome and exposure se)

    do_single_term_mr_estimate(self, exposure_tuple, outcome_tuple)
        single term MR estimate. In contrast to SMR estimate, uses a single variance term.

    do_single_term_mr_estimate(self, exposure_tuple, outcome_tuple)
        single term MR estimate. In contrast to SMR estimate, uses a single variance term.
        adds it to the list of estimates.

    mr_presso(self, n_sims=1000, significance_thresh=0.05):
        Does MR-PRESSO, number of sims is the number of permutations done, and the significance term is the term
        at which estimates are rejected.

    do_lda_mr_egger(self, ld_matrix)
        Does LDA-MR-Egger.
        requires an LD matrix (peason correlation), which is ordered (rows and columns) by the estimates that were made.

    do_lda_mr_egger_on_estimates(self, list_of_outcome_tuples, list_of_exposure_tuples, pearson_ld_matrix,
                                write_out=False):
        Does LDA MR Egger
        Requires the list of (beta, se) tuples from the outcome and exposure and the pearson_ld_matrix.
        write_out is for debug purposes.


    """

    def __init__(self):

        self.estimation_done = False  # when the estimation is called, this will return true
        self.estimation_data = []     # will be filled with tuples containing betas and standard error
        self.estimation_snps = []

        # these are optional

        self.outcome_tuples = []
        self.exposure_tuples = []

        #cochrans q test
        self.q_test_indices_remaining = None


        # this will be used to do meta analysis techniques.
        # meaning there will be some
        self._ivw_intermediate_top = []
        self._ivw_intermediate_bottom = []




    def add_estimate(self, beta_se_tuple, variant_name, pos, chr):
        """
        Adds an estimate to the class.

        :param beta_se_tuple: Estimation data.
        :param variant_name: name of variant for external reference
        :param pos: position of the variant for external reference
        :param chr: chromosome of the variant for external reference
        :return:
        """
        self.estimation_done = False

        self.estimation_data.append(beta_se_tuple)
        self.estimation_snps.append((variant_name, pos, chr))

    def do_ivw_heterogeneity_estimation(self):
        """
        Run a heterogeneity estimation on an already done IVW estimation.
        :return: a p value of estimation.
        """
        heterogeneity_p_value = self.do_ivw_heterogeneity_estimation_on_estimate_vector(self.estimation_data)
        return heterogeneity_p_value

    def do_ivw_heterogeneity_estimation_on_estimate_vector(self, estimation_vector):
        """
        Provides a heterogeneity estimation p value of an ivw estimation based on an estimation vector.

        :param estimation_vector: estimation vector of MR estimates.
        :return: a p value of heterogeneity.
        """
        if len(self.estimation_data) < 3:
            raise ValueError("Less than three estimates supplied, cannot do cochrans q heterogeneity analysis")

        ivw_estimate = self.do_ivw_estimation_on_estimate_vector(estimation_vector, save_intermediate=False)

        chi_sq_sum = 0.0
        for sub_estimate in estimation_vector:
            chi_sq_sum += (ivw_estimate[0] - sub_estimate[0]) **2 / sub_estimate[1] ** 2

        return scipy.stats.chi2.sf(chi_sq_sum, len(estimation_vector) - 1)


    def do_ivw_estimation(self):
        """
        Does IVW estimation on all the methods

        :return: tuple of floats: beta, se, wald_p_val of the estimate.
        """

        beta, se, wald_p_val = self.do_ivw_estimation_on_estimate_vector(
            self.estimation_data, save_intermediate=True)
        self.estimation_done = True

        return beta, se, wald_p_val


    def do_ivw_estimation_on_estimate_vector(self, estimation_vec, save_intermediate=False):
        """
        Estimates IVW on a specified estimate vector

        :param estimation_vec: list of (beta,se) tuples
        :param save_intermediate: save intermediate results to the class.
        :return: tuple of floats: beta, se, wald_p_val of the estimate.
        """
        if len(estimation_vec) == 0:
            raise RuntimeError('No estimates supplied to do estimation')

        beta_top = 0.0
        beta_bottom = 0.0

        ivw_intermediate_top = []
        ivw_intermediate_bottom = []

        i = 0

        for smr_result in estimation_vec:

            # make sure the standard error cannot be zero.
            if smr_result[1] == 0:
                smr_result[1] = np.nextafter(0, 1)


            ivw_intermediate_top.append(smr_result[0] * (smr_result[1] ** -2))
            ivw_intermediate_bottom.append((smr_result[1] ** -2))



            beta_top += ivw_intermediate_top[i]
            beta_bottom += ivw_intermediate_bottom[i]
            i += 1

        if save_intermediate:
            self._ivw_intermediate_top = ivw_intermediate_top
            self._ivw_intermediate_bottom = ivw_intermediate_bottom

        beta_ivw = beta_top / beta_bottom
        se_ivw = math.sqrt(1 / beta_bottom)

        p_value = scipy.stats.norm.sf(abs(beta_ivw / se_ivw)) * 2

        return beta_ivw, se_ivw, p_value

    def get_ivw_estimates(self):
        """
        Mirrors do_ivw_estimation

        :return: tuple of floats: beta, se, wald_p_val of the estimate.
        """
        warnings.warn("This method is deprecated, please use the do_ivw_estimation() method.", DeprecationWarning)
        return self.do_ivw_estimation()

    def do_smr_estimate(self, exposure_tuple, outcome_tuple):
        """
        Determine SMR test effect and standard error.
        Identifies standard error of the estimate, based on on two variance terms

        :param exposure_data: Exposure estimates which must have the methods get_beta and get_z_score
        :param outcome_data: Outcome estimates which must have the methods get_beta and get_z_score
        :return: tuple of smr beta and smr se of the estimate.
        """

        z_score_exposure = exposure_tuple[0] / exposure_tuple[1]
        z_score_outcome = outcome_tuple[0] / outcome_tuple[1]


        # from SMR paper.
        t_stat = ((z_score_exposure ** 2) * (z_score_outcome ** 2)) \
                 / \
                 ((z_score_exposure ** 2) + (z_score_outcome ** 2))

        # this needs a check,
        # checked it with results from Zhy et al. see validate_SMR_estimates.py
        p_value = scipy.stats.chi2.sf(t_stat, 1)
        z_score = scipy.stats.norm.isf(p_value / 2)

        # also from Zhu et al.
        beta_smr = outcome_tuple[0] / exposure_tuple[0]
        se_smr = abs(beta_smr / z_score)

        return [beta_smr, se_smr, p_value]

    def do_and_add_smr_estimation(self, exposure_tuple, outcome_tuple, variant_name=None, pos=None, chr=None):
        """
        Does an SMR estimation (two variance terms included) and adds it to the class.


        :param exposure_tuple: beta,se tuple of the exposure summary statistics
        :param outcome_tuple: beta, se tuple of the outcome summary statistics
        :param variant_name: name of the variant for external reference
        :param pos: position of the variant for external reference
        :param chr: chromosome of the variant for external reference
        :return:
        """


        estimates = self.do_smr_estimate(exposure_tuple, outcome_tuple)
        self.estimation_data.append(estimates)
        self.estimation_snps.append((variant_name, pos, chr))

        self.outcome_tuples.append(outcome_tuple)
        self.exposure_tuples.append(exposure_tuple)

    def do_weighted_median_meta_analysis(self):
        return self.do_weighted_median_meta_analysis_on_estimate_vectors(self.exposure_tuples, self.outcome_tuples)

    def do_weighted_median_meta_analysis_on_estimate_vectors(self, exposure_associations, outcome_associations):
        """
        Performs the weighted median estimation
        following Web Appendix 2 from
        Genet Epidemiol. 2016 May; 40(4): 304–314.
        Published online 2016 Apr 7. doi: 10.1002/gepi.21965
        PMCID: PMC4849733
        PMID: 27061298
        Consistent Estimation in Mendelian Randomization with Some Invalid Instruments Using a Weighted Median Estimator
        Jack Bowden, 1 George Davey Smith, 1 Philip C. Haycock, 1 and Stephen Burgess

        :param estimation_vector:
        :return: Weighted median estimator.

        """

        def weighted_median_beta(betas_per_iv, weights):
            if len(betas_per_iv) < 2:
                raise ValueError("Weighted median estimator requires a minimum of 2 ivs")
            if len(weights) != len(betas_per_iv):
                raise ValueError("Betas and weights should have the same shape")

            ordering = np.argsort(betas_per_iv)
            betas = betas_per_iv[ordering]
            weights = 1 / (weights[ordering] ** 2)
            weight_sum = (np.cumsum(weights) - 0.5* weights) / np.sum(weights)
            below = np.max(np.where(weight_sum < 0.5)[0])
            weighted_estimate = betas[below] + (betas[below+1]-betas[below]) * (
                                0.5-weight_sum[below]) / (weight_sum[below+1] - weight_sum[below])

            return weighted_estimate

        def bootstrap_se_of_wm_estimate(exposure_tuple, outcome_tuple, weights, n_bootstraps=1000):

            exposure_mat = np.asarray(exposure_tuple)[:,:2]
            outcome_mat = np.asarray(outcome_tuple)[:,:2]

            exposure_draws = np.random.normal(exposure_mat[:,0], exposure_mat[:,1],
                                              size=(n_bootstraps, exposure_mat.shape[0]))

            outcome_draws = np.random.normal(outcome_mat[:,0], outcome_mat[:,1],
                                             size=(n_bootstraps, outcome_mat.shape[0]))

            mr_draws = outcome_draws / exposure_draws

            bootstrapped_wm = np.zeros(n_bootstraps)
            for i in range(n_bootstraps):
                bootstrapped_wm[i] = weighted_median_beta(mr_draws[i,:], weights)
            return np.std(bootstrapped_wm)

        if len(exposure_associations) != len(outcome_associations):
            raise ValueError(f"Exposure tuples should have the same length as the outcome tuple")

        estimates  = [self.do_single_term_mr_estimate(exposure_associations[i], outcome_associations[i]) for i in
                              range(len(exposure_associations))]

        mr_betas = np.asarray([x[0] for x in estimates])
        mr_weights = np.asarray([x[1] if x[1] > 1e-10 else 1e-10 for x in estimates])

        beta_wm = weighted_median_beta(mr_betas, mr_weights)
        se_wm = bootstrap_se_of_wm_estimate(exposure_associations, outcome_associations, mr_weights)
        p_val = scipy.stats.t.sf(abs(beta_wm / se_wm), df = len(mr_betas) - 1)

        return beta_wm, se_wm, p_val


    def do_chochrans_q_meta_analysis(self, p_value_threshold):
        """
        Does a chochrans Q meta analysis on the interally present estimates, and estimates a combined causal effect.

        :param p_value_threshold: p value threshold when to reject the null hypothesis that the estimate is drawn from
        the same distribution.

        :return: tuple of floats: beta, se, wald_p_val of the estimate after chochran's Q meta analysis.
        """
        if len(self.estimation_data) < 2:
            raise ValueError("Less than 2 estimates supplied, cannot do cochrans q analysis")

        if not self.estimation_done:
            raise ValueError("No causal_inference estimation done.")


        #base the calculations on indices that are remaining
        indices_remaining = list(range(len(self.estimation_data)))
        old_indices = indices_remaining

        #these will change while the algorithm is running.
        beta, se, wald_p_val = self.do_ivw_estimation()

        tmp_beta_ivw = beta
        tmp_se_ivw = se
        tmp_p_ivw = wald_p_val
        save_beta_ivw, save_se_ivw, save_p_ivw = tmp_beta_ivw, tmp_se_ivw, tmp_p_ivw

        p_val = 0.0

        while (p_value_threshold > p_val) and (len(indices_remaining) >= 2):
            old_indices = copy.deepcopy(indices_remaining)

            #determine the Q statistic, per estimate.
            q_terms = [
                self._ivw_intermediate_bottom[i] * (self.estimation_data[i][0] - tmp_beta_ivw) ** 2
                for i in indices_remaining
            ]

            # determine the total Q statistic and find the p value.
            q_stat = np.sum(q_terms)

            max_stat = max(q_terms)
            max_indice = q_terms.index(max_stat)

            p_val = scipy.stats.chi2.sf(q_stat,  len(q_terms)-1)


            del indices_remaining[max_indice]

            save_beta_ivw, save_se_ivw, save_p_ivw = tmp_beta_ivw, tmp_se_ivw, tmp_p_ivw

            #determine the new causal_inference estimates
            tmp_beta_ivw, tmp_se_ivw, tmp_p_ivw = \
                self.do_ivw_estimation_on_estimate_vector([self.estimation_data[i] for i in indices_remaining])

        #save it now.
        self.q_test_indices_remaining = old_indices


        return (save_beta_ivw, save_se_ivw, save_p_ivw), old_indices, p_val

    def do_egger_regression(self):
        """
        Does egger regression based on single variance term estimates.

        :return: list of length two each with a tuple of floats: beta, se, wald_p_val of the estimate for intercept
        and slope
        """
        return self.do_egger_regression_single_variance_term()

    def do_egger_regression_single_variance_term(self):
        """
        Does egger regression based on single variance term estimates.

        :return: list of length two each with a tuple of floats: beta, se, wald_p_val of the estimate for intercept
        and slope respectively.
        """

        num_estimates = len(self.estimation_data)

        # runtime checks.
        if num_estimates < 3:
            raise ValueError("Only {} estimates supplied, need at least three to estimate egger".format(num_estimates))

        if len(self.exposure_tuples) != num_estimates:
            raise ValueError("No exposure data present, cannot do Egger regression.")

        if len(self.outcome_tuples) != num_estimates:
            raise ValueError("No outcome data present, cannot do Egger regression.")

        """
        Now turn exposure into positive values.
        """
        outcome_tuples = copy.deepcopy(self.outcome_tuples)
        exposure_tuples = copy.deepcopy(self.exposure_tuples)

        for i in range(num_estimates):
            if exposure_tuples[i][0] < 0:
                # flip.
                exposure_tuples[i] = (-1.0 * exposure_tuples[i][0], exposure_tuples[i][0])
                outcome_tuples[i] = (-1.0 * outcome_tuples[i][0], outcome_tuples[i][0])

        x_dat = np.asarray([x[0] for x in exposure_tuples], dtype=float)
        x_dat = add_constant(x_dat)

        y_dat = np.asarray([x[0] for x in outcome_tuples], dtype=float)

        w_dat = np.zeros(len(self.estimation_data), dtype=float)

        for i in range(len(self.estimation_data)):
            w_dat[i] = 1 / (float(outcome_tuples[i][1]) ** 2)

        wls_model = WLS(y_dat, x_dat, weights=w_dat)
        results = wls_model.fit()


        self.egger_intercept = (results.params[0], results.bse[0], results.pvalues[0])
        self.egger_slope = (results.params[1], results.bse[1], results.pvalues[1])

        self.egger_done = True

        return self.egger_intercept, self.egger_slope


    def do_egger_regression_two_variance_term(self):

        """
        Does egger regression based on two variance term estimates.

        :return: list of length two each with a tuple of floats: beta, se, wald_p_val of the estimate for intercept
        and slope respectively.
        """


        num_estimates = len(self.estimation_data)

        # runtime checks.

        if num_estimates < 3:
            raise ValueError("Only {} estimates supplied, need at least three to estimate egger".format(num_estimates))

        if len(self.exposure_tuples) != num_estimates:
            raise ValueError("No exposure data present, cannot do Egger regression.")

        if len(self.outcome_tuples) != num_estimates:
            raise ValueError("No outcome data present, cannot do Egger regression.")

        """
        Now turn exposure into positive values.
        """

        outcome_tuples = copy.deepcopy(self.outcome_tuples)
        exposure_tuples = copy.deepcopy(self.exposure_tuples)

        for i in range(num_estimates):
            if exposure_tuples[i][0] < 0:
                # flip.
                exposure_tuples[i] = (-1*exposure_tuples[i][0],exposure_tuples[i][0])
                outcome_tuples[i] = (-1*outcome_tuples[i][0], outcome_tuples[i][0])


        x_dat = np.asarray([x[0] for x in exposure_tuples])
        x_dat = add_constant(x_dat)

        y_dat = np.asarray([x[0] for x in  outcome_tuples])

        #if this value is zero, we add the smallest possible constant, so it can still be used as weights.
        #checked with the 2015 paper introducing MR-Egger, and it works as expected.

        w_dat = np.zeros(len(self.estimation_data))
        for i in range(len(self.estimation_data)):
            w_dat[i] = outcome_tuples[i][0] ** -2  / \
                       ( (outcome_tuples[i][0]**-2 * outcome_tuples[i][1] ** 2) +
                         (exposure_tuples[i][0]**-2 * exposure_tuples[i][1] ** 2)
                    )


        wls_model = WLS(y_dat, x_dat, weights=w_dat)
        results = wls_model.fit()

        self.egger_intercept = (results.params[0], results.bse[0], results.pvalues[0])
        self.egger_slope = (results.params[1], results.bse[1], results.pvalues[1])

        self.egger_done = True

        return self.egger_intercept, self.egger_slope


    def do_single_term_mr_estimate(self, exposure_tuple, outcome_tuple):
        """
        Does a single variance term MR estimate on a variant.

        :param exposure_tuple: beta, se tuple of the exposure
        :param outcome_tuple: beta, se tuple of the outcome
        :return: beta, se and p value of the single estimate.
        """
        beta_mr = outcome_tuple[0] / exposure_tuple[0]
        se_mr = np.sqrt((outcome_tuple[1] ** 2) / (exposure_tuple[0] ** 2))

        z_score = beta_mr / se_mr
        p_value = scipy.stats.norm.sf(abs(z_score)) * 2

        return beta_mr, se_mr, p_value


    def do_and_add_single_term_mr_estimation(self, exposure_tuple, outcome_tuple):
        """
        Does a single term variance estimate of a variant, and adds it to the class for further analysis.

        :param exposure_tuple: beta, se tuple of the exposure
        :param outcome_tuple: beta, se tuple of the outcome
        :return: None
        """

        estimates = self.do_single_term_mr_estimate(exposure_tuple, outcome_tuple)
        self.estimation_data.append(estimates)

        self.outcome_tuples.append(outcome_tuple)
        self.exposure_tuples.append(exposure_tuple)

    def mr_presso(self, n_sims=1000, significance_thresh=0.05):
        """
        Python reimplementation of MR-PRESSO.

        :param n_sims: number of permutation simulations.
        :param significance_thresh: significance thresshold.
        :return: beta, se and p value of the estimate after the bad snps were removed. If no estimate can be made,
        returns a tuple of 3* (np.nan).
        """

        def make_random_data():
            beta_ivw, _, _ = self.do_ivw_estimation()
            random_exposure = np.random.normal([x[0] for x in self.exposure_tuples],
                                               [x[1] for x in self.exposure_tuples])
            random_outcome = np.random.normal([beta_ivw * x[0] for x in self.exposure_tuples],
                                              [x[1] for x in self.outcome_tuples])

            mr_estimates = np.zeros((len(random_outcome), 3), dtype=float)
            for i in range(len(random_outcome)):
                mr_estimates[i, :] = self.do_single_term_mr_estimate(
                    (random_exposure[i], self.exposure_tuples[i][1]),
                    (random_outcome[i], self.outcome_tuples[i][1]))

            return random_exposure, random_outcome, mr_estimates

        def leave_one_out_residual_sum_of_squares(estimation_data,
                                                  weighted_outcome,
                                                  weighted_exposure):

            estimation_data = np.asarray(estimation_data)
            leave_one_out_ivw = np.zeros(shape=(len(estimation_data), 3), dtype=float)
            for i in range(len(estimation_data)):
                leave_one_out_ivw[i, :] = self.do_ivw_estimation_on_estimate_vector(
                    np.delete(estimation_data, i, 0)
                )

            rss = (weighted_outcome - leave_one_out_ivw[:, 0] * weighted_exposure) ** 2

            return rss, leave_one_out_ivw

        def make_random_data_and_return_rss(weights):
            exposure, outcome, mr_estimates = make_random_data()

            weighted_exposure = exposure * weights
            weighted_outcome = outcome * weights

            rss, _ = leave_one_out_residual_sum_of_squares(mr_estimates, weighted_outcome, weighted_exposure)

            return np.sum(rss), np.concatenate(
                (exposure.reshape(len(exposure), 1), outcome.reshape(len(outcome), 1)), axis=1)

        def randomly_sample_distortion(outlier_indices):
            estimates = np.asarray(self.estimation_data)
            estimates_no_outliers = np.delete(estimates, outlier_indices, axis=0)
            estimates_only_outliers = estimates[outlier_indices, :][0]

            indices_sampled_from_no_outliers = np.random.choice(estimates_no_outliers.shape[0],
                                                                size=estimates_no_outliers.shape[0],
                                                                replace=True)
            return self.do_ivw_estimation_on_estimate_vector(
                np.concatenate(
                    (estimates_no_outliers[indices_sampled_from_no_outliers, :], estimates_only_outliers))
            )

        # runtime checks.
        num_estimates = len(self.estimation_data)

        if num_estimates < 3:
            raise ValueError(
                "Only {} estimates supplied, need at least three to find simulate_mr presso outliers".format(
                    num_estimates))

        if len(self.exposure_tuples) != num_estimates:
            raise ValueError("No exposure sumstats present, cannot do mr_presso outlier.")

        if len(self.outcome_tuples) != num_estimates:
            raise ValueError("No outcome sumstats present, cannot do mr_presso outlier.")

        # this is just following MR presso.
        outcome = np.asarray(self.outcome_tuples, dtype=float)
        exposure = np.asarray(self.exposure_tuples, dtype=float)
        weighted_outcome = np.asarray([x[0] / np.sqrt(x[1] ** 2) for x in self.outcome_tuples], dtype=float)
        weighted_exposure = np.asarray([self.exposure_tuples[i][0] / np.sqrt(self.outcome_tuples[i][1] ** 2)
                                        for i in range(len(self.exposure_tuples))], dtype=float)
        weights = np.asarray([1 / np.sqrt(x[1] ** 2) for x in self.outcome_tuples], dtype=float)

        rss, list_of_assocs = leave_one_out_residual_sum_of_squares(self.estimation_data,
                                                                    weighted_outcome,
                                                                    weighted_exposure)

        expected_results = [make_random_data_and_return_rss(weights) for _ in range(n_sims)]

        sim_rss = [x[0] for x in expected_results]

        global_p_val = sum(sim_rss > sum(rss)) / n_sims
        local_p_val = None
        if global_p_val < significance_thresh:
            expected_betas = np.zeros((num_estimates, n_sims, 2), dtype=float)
            for i in range(n_sims):
                expected_betas[:, i] = expected_results[i][1]

            difference = outcome[:, 0] - exposure[:, 0] * list_of_assocs[:, 0]
            expected_difference = expected_betas[:, :, 1] - expected_betas[:, :, 0] * np.tile(list_of_assocs[:, 0],
                                                                                              (n_sims,
                                                                                               1)).transpose()
            local_p_val = np.sum(expected_difference ** 2 > (difference ** 2).reshape((len(difference), 1)),
                                 axis=1) / n_sims
            local_p_val = np.asarray(
                [x * len(difference) if x * len(difference) < 1.0 else 1.0 for x in local_p_val])

        # distortion test.
        outlier_corrected_ivw_result = (np.nan, np.nan, np.nan)
        if local_p_val is not None and sum(local_p_val < significance_thresh):
            outliers = local_p_val < significance_thresh

            exposure_betas = [self.exposure_tuples[i][0] for i in range(len(self.estimation_data)) if not outliers[i]]
            outcome_betas = [self.outcome_tuples[i][0] for i in range(len(self.estimation_data)) if not outliers[i]]
            weights = [1 / self.outcome_tuples[i][1] ** 2 for i in range(len(self.estimation_data)) if not outliers[i]]

            outlier_corrected_ivw_result = WLS(exog=exposure_betas, endog=outcome_betas, weights=weights).fit()

            return outlier_corrected_ivw_result.params[0], outlier_corrected_ivw_result.bse[0], outlier_corrected_ivw_result.pvalues[0]
        else:
            return outlier_corrected_ivw_result


    def do_lda_mr_egger(self, ld_matrix):
        """
        Perform LDA MR egger on internal estimates.
        :param ld_matrix: pearson LD matrix. ordered by the estimates.
        :return: list of length two each with a tuple of floats: beta, se, wald_p_val of the estimate for intercept
        and slope respectively.
        """
        return self.do_lda_mr_egger_on_estimates(self.outcome_tuples, self.exposure_tuples, ld_matrix)



    def do_lda_mr_egger_on_estimates(self, list_of_outcome_tuples, list_of_exposure_tuples, pearson_ld_matrix, write_out=False):
        """
        This will do LDA simulate_mr egger regression as described in Barfield et al. 2018, genetic epidemiology
        Implemented based on their paper, and a reference implementation they provided in personal communication

        <Begin email.>
        Hi Adriaan,
        Below please find the R function to implement the approach. Please let me know if you have any questions.

        -Richard

        X is the vector of  joint eQTL effects
        Y is the vector of joint GWAS effects
        W is the inverse of the covariance of the joint GWAS effects (i.e. var(Y))

        weight.func2<-function(X,Y,W){
          bX<-cbind(1,X)
          bread<-solve(crossprod(bX,W)%*%bX)
          theEsts<-bread%*%crossprod(bX,W%*%Y)
          theresid<-c(Y-theEsts[1]-X*theEsts[2])
          Sig.Est<-c(crossprod(theresid,W%*%theresid))/(length(X)-2)
          finresults<- cbind(theEsts,diag(bread)*Sig.Est)
          TestStat<-theEsts/sqrt(finresults[,2])
          Pvals<-2*pt(abs(TestStat),df = nrow(bX)-2,lower.tail = F)
          return(cbind(finresults,TestStat,Pvals))
        }
        <End email.>


        :param list_of_outcome_tuples: list of beta,se tuples from the outcome
        :param list_of_exposure_tuples: list of beta,se tuples from the exposure
        :param pearson_ld_matrix: peason LD matrix in the order of the estimates.
        :return: :return: list of length two each with a tuple of floats: beta, se, wald_p_val of the estimate for intercept
        and slope respectively.
        """


        if len(list_of_outcome_tuples) < 3:
            raise ValueError("Could not do lda simulate_mr egger on estimates, too little estimates supplied")

        marginal_exposure = np.asarray(list_of_exposure_tuples, dtype=float)
        marginal_outcome = np.asarray(list_of_outcome_tuples, dtype=float)

        #flip to make exposure strictly positive.
        to_flip = marginal_exposure[:, 0] < 0.0

        marginal_exposure[to_flip, 0] = marginal_exposure[to_flip,0] * -1
        marginal_outcome[to_flip, 0] = marginal_outcome[to_flip, 0] * -1
        for i in to_flip:
            if i:
                pearson_ld_matrix[i, :] = pearson_ld_matrix[i, :] * -1
                pearson_ld_matrix[:, i] = pearson_ld_matrix[:, i] * -1

        if marginal_exposure.shape[1] < 1:
            raise ValueError("No standard errors supplied to the marginal exposure")

        if marginal_outcome.shape[1] < 1:
            raise ValueError("No standard errors supplied to the marginal outcome")

        sigma = pearson_ld_matrix

        inv_sigma = np.linalg.inv(sigma)

        conditional_outcome = inv_sigma @ marginal_outcome[:,0]

        conditional_exposure = inv_sigma @ marginal_exposure[:,0]

        sigma_g = inv_sigma

        b_x = np.concatenate((np.ones((conditional_exposure.shape[0],1)) , conditional_exposure.reshape(conditional_exposure.shape[0], 1)), axis=1)

        bread = np.linalg.inv(b_x.transpose() @ sigma_g @ b_x)

        estimates = bread @ b_x.transpose() @ sigma_g @ conditional_outcome

        residuals = conditional_outcome - estimates[0] - conditional_exposure * estimates[1]

        significant_estimates = (residuals.transpose() @ (sigma_g @ residuals)) / (conditional_exposure.shape[0] - 2)

        test_stat = estimates / np.sqrt(significant_estimates * np.diag(bread))

        p_val = 2 * scipy.stats.t.sf(np.abs(test_stat), df=conditional_exposure.shape[0]-2)

        """
        This only used to compare to the R implementation.
        """
        if write_out:
            with open("check_r_implementation.txt", "w") as f:
                for i in range(marginal_outcome.shape[0]):
                    f.write("{}\t{}\t".format(conditional_exposure[i], conditional_outcome[i]) + "\t".join([str(x) for x in sigma_g[i,:]]) + "\n" )

        return estimates, np.sqrt(np.diag(bread) * significant_estimates), test_stat, p_val