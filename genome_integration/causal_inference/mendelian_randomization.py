"""
This class is used to do a 'simple'
inverse variance estimate on independent effects.
"""

__author__      = "Adriaan van der Graaf"
__copyright__   = "Copyright 2018, Adriaan van der Graaf"

import math
import copy
import numpy as np
import scipy.stats
from statsmodels.regression.linear_model import WLS
from statsmodels.tools.tools import add_constant


from .. import utils
from .. import association




class MendelianRandomization(association.BaseAssociation):

    def __init__(self):
        super().__init__()

        self.estimation_done = False  # when the estimation is called, this will return true
        self.estimation_data = []     # will be filled with tuples containing betas and standard error
        self.estimation_snps = []

        # these are optional

        self.outcome_tuples = []
        self.exposure_tuples = []

        # this will be used to do meta analysis techniques.
        # meaning there will be some
        self.ivw_intermediate_top = []
        self.ivw_intermediate_bottom = []

        self.beta = np.nan
        self.se = np.nan
        self.wald_p_val = np.nan

        #cochrans q done.
        self.q_test_done = False
        self.q_test_p_value = None
        self.q_test_indices_remaining = None

        self.q_test_beta = None
        self.q_test_se = None
        self.q_test_p_value = None

        #egger regression
        self.egger_done = False

        self.egger_intercept = None
        self.egger_slope = None



    def write_for_smr_style_plot_without_q_data(self, filename):
        """
        will write the cross style plot.


        :param filename:
        :return:
        """
        with open(filename, "w") as f:
            f.write("snp_name\tbeta_exposure\tse_exposure\tbeta_outcome\tse_outcome\tbeta_smr\tse_ivw\n")
            for i in range(len(self.estimation_data)):
                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.estimation_snps[i][0],
                                                              self.exposure_tuples[i][0],
                                                              self.exposure_tuples[i][1],
                                                              self.outcome_tuples[i][0],
                                                              self.outcome_tuples[i][1],
                                                              self.estimation_data[i][0],
                                                              self.estimation_data[i][1]))



    def write_for_smr_style_plot_with_q_data(self, filename):
        """
        Contains an extra column.

        :param filename:
        :return:
        """
        with open(filename, "w") as f:
            f.write("snp_name\tbeta_exposure\tse_exposure\tbeta_outcome\tse_outcome\tbeta_smr\tse_ivw\tq_pruned\n")
            for i in range(len(self.estimation_data)):
                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.estimation_snps[i][0],
                                                              self.exposure_tuples[i][0],
                                                              self.exposure_tuples[i][1],
                                                              self.outcome_tuples[i][0],
                                                              self.outcome_tuples[i][1],
                                                              self.estimation_data[i][0],
                                                              self.estimation_data[i][1],
                                                              0.6 + 0.4 * float(i in self.q_test_indices_remaining),
                                                              ))


    def write_for_smr_style_plot(self, filename):

        if self.q_test_done:
            self.write_for_smr_style_plot_with_q_data(filename)
        else:
            self.write_for_smr_style_plot_without_q_data(filename)



    def add_estimate(self, beta_se_tuple, snp_name_1, pos, chr, snp_name_2):
        self.estimation_done = False
        self.estimation_data.append(beta_se_tuple)
        self.estimation_snps.append((snp_name_1, pos, chr, snp_name_2))

    def do_ivw_estimation(self):

        self.beta, self.se, self.wald_p_val = self.do_ivw_estimation_on_estimate_vector(
            self.estimation_data, self.ivw_intermediate_top, self.ivw_intermediate_bottom)
        self.estimation_done = True
        return self.beta, self.se, self.wald_p_val

    def do_ivw_estimation_on_estimate_vector(self, estimation_vec, ivw_intermediate_top=None, ivw_intermediate_bottom=None):

        if len(estimation_vec) == 0:
            raise RuntimeError('No estimates supplied to do estimation')

        beta_top = 0.0
        beta_bottom = 0.0

        if ivw_intermediate_top is None or ivw_intermediate_bottom is None:
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

        beta_ivw = beta_top / beta_bottom
        se_ivw = math.sqrt(1 / beta_bottom)

        p_value = scipy.stats.norm.sf(abs(beta_ivw / se_ivw)) * 2

        return beta_ivw, se_ivw, p_value

    def get_ivw_estimates(self):
        if self.estimation_done:
            return self.beta, self.se, self.wald_p_val
        else:
            return self.do_ivw_estimation()

    def do_smr_estimate(self, exposure_tuple, outcome_tuple):
        """
        Determine SMR test effect and standard error.

        :param exposure_data: Exposure estimates which must have the methods get_beta and get_z_score
        :param outcome_data: Outcome estimates which must have the methods get_beta and get_z_score
        :return: tuple of smr beta and smr se
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

    def do_and_add_smr_estimation(self, exposure_tuple, outcome_tuple, snp_name_1, pos, chr, snp_name_2):
        estimates = self.do_smr_estimate(exposure_tuple, outcome_tuple)
        self.estimation_data.append(estimates)
        self.estimation_snps.append((snp_name_1, pos, chr, snp_name_2))

        self.outcome_tuples.append(outcome_tuple)
        self.exposure_tuples.append(exposure_tuple)


    def write_ivw_header(self):
        return "exposure\toutcome\tchromosome\tn_snps\tbeta_ivw\tse_ivw\tp_ivw\testimation_snp_names\testimation_snp_names\testimation_snp_betas\testimation_snp_se"

    def write_ivw_line(self, exposure, outcome, chromosome):
        if not self.estimation_done:
            return None

        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            exposure,
            outcome,
            chromosome, #this requires that the estimation snps are well defined.
            len(self.estimation_data),
            self.beta,
            self.se,
            self.wald_p_val,
            ','.join([x[0] for x in self.estimation_snps]),
            ','.join([str(x[0]) for x in self.estimation_data]),
            ','.join([str(x[1]) for x in self.estimation_data])
        )

    def write_ivw_estimates(self, filename):
        strings = ["beta_ivw\tse_ivw\tp_val_ivw\testimate_type\tsnp_1\tbp\tchr\tsnp_2",
                   "\t".join(
                            [
                            str(self.beta),
                            str(self.se),
                            str(self.wald_p_val),
                            "ivw_estimate",
                            "NA",
                            "NA",
                            "NA",
                            "NA"
                            ]
                            )
                   ]

        for i in range(len(self.estimation_data)):
            strings.append(
                "\t".join(
                        [str(self.estimation_data[i][0]),
                         str(self.estimation_data[i][1]),
                         str(self.estimation_data[i][2]),
                         "smr_estimate",
                         str(self.estimation_snps[i][0]),
                         str(self.estimation_snps[i][1]),
                         str(self.estimation_snps[i][2]),
                         str(self.estimation_snps[i][3])
                        ]
                        )
            )

        utils.write_list_to_newline_separated_file(strings, filename)


    def do_chochrans_q_meta_analysis(self, p_value_threshold):

        if len(self.estimation_data) < 3:
            raise ValueError("Less than three estimates supplied, cannot do cochrans q analysis")

        if not self.estimation_done:
            raise ValueError("No causal_inference estimation done.")


        #base the calculations on indices that are remaining
        indices_remaining = list(range(len(self.estimation_data)))
        old_indices = indices_remaining

        #these will change while the algorithm is running.
        tmp_beta_ivw = self.beta
        tmp_se_ivw = self.se
        tmp_p_ivw = self.wald_p_val
        save_beta_ivw, save_se_ivw, save_p_ivw = tmp_beta_ivw, tmp_se_ivw, tmp_p_ivw

        p_val = 0.0

        while (p_value_threshold > p_val) and (len(indices_remaining) >= 3):
            old_indices = copy.deepcopy(indices_remaining)

            #determine the Q statistic, per estimate.
            q_terms = [
                self.ivw_intermediate_bottom[i] * (self.estimation_data[i][0] - tmp_beta_ivw) ** 2
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
        self.q_test_done = True
        self.q_test_p_value = p_value_threshold
        self.q_test_indices_remaining = old_indices

        self.q_test_beta = save_beta_ivw
        self.q_test_se = save_se_ivw
        self.q_test_p_value = save_p_ivw

        return (save_beta_ivw, save_se_ivw, save_p_ivw), old_indices, p_val

    def do_egger_regression(self):
        return self.do_egger_regression_single_variance_term()

    def do_egger_regression_single_variance_term(self):

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
                exposure_tuples[i] = (-1 * exposure_tuples[i][0], exposure_tuples[i][0])
                outcome_tuples[i] = (-1 * outcome_tuples[i][0], outcome_tuples[i][0])

        x_dat = np.asarray([x[0] for x in exposure_tuples])
        x_dat = add_constant(x_dat)

        y_dat = np.asarray([x[0] for x in outcome_tuples])


        w_dat = np.zeros(len(self.estimation_data))
        for i in range(len(self.estimation_data)):
            w_dat[i] = 1 / (outcome_tuples[i][1] ** 2)

        wls_model = WLS(y_dat, x_dat, weights=w_dat)
        results = wls_model.fit()

        self.egger_intercept = (
        results.params[0], results.bse[0], scipy.stats.norm.sf(abs(results.params[0] / results.bse[0])))
        self.egger_slope = (
        results.params[1], results.bse[1], scipy.stats.norm.sf(abs(results.params[1] / results.bse[1])))

        self.egger_done = True

        return self.egger_intercept, self.egger_slope




    def do_egger_regression_two_variance_term(self):

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
                       ( (outcome_tuples[i][0]**-2 * outcome_tuples[i][1] **2) +
                         (exposure_tuples[i][0]**-2 * exposure_tuples[i][1] **2)
                    )


        wls_model = WLS(y_dat, x_dat, weights=w_dat)
        results = wls_model.fit()

        self.egger_intercept = (results.params[0], results.bse[0], scipy.stats.norm.sf(abs(results.params[0]/results.bse[0])))
        self.egger_slope = (results.params[1], results.bse[1], scipy.stats.norm.sf(abs(results.params[1]/results.bse[1])))

        self.egger_done = True

        return self.egger_intercept, self.egger_slope


    def do_single_term_mr_estimate(self, exposure_tuple, outcome_tuple):
        """
        Determine SMR test effect and standard error.

        :param exposure_data: Exposure estimates which must have the methods get_beta and get_z_score
        :param outcome_data: Outcome estimates which must have the methods get_beta and get_z_score
        :return: tuple of smr beta and smr se
        """

        beta_mr = outcome_tuple[0] / exposure_tuple[0]
        se_mr = np.sqrt((outcome_tuple[1] ** 2) / (exposure_tuple[0] ** 2))

        z_score = beta_mr / se_mr
        p_value = scipy.stats.norm.sf(abs(z_score)) * 2

        return [beta_mr, se_mr, p_value]

    def do_and_add_single_term_mr_estimation(self, exposure_tuple, outcome_tuple):
        estimates = self.do_single_term_mr_estimate(exposure_tuple, outcome_tuple)
        self.estimation_data.append(estimates)

        self.outcome_tuples.append(outcome_tuple)
        self.exposure_tuples.append(exposure_tuple)

    def mr_presso(self, n_sims=1000, significance_thresh=0.05):

        def make_random_data():
            beta_ivw, _, _ = self.get_ivw_estimates()
            random_exposure = np.random.normal([x[0] for x in self.exposure_tuples],
                                               [x[1] for x in self.exposure_tuples])
            random_outcome = np.random.normal([beta_ivw * x[0] for x in self.exposure_tuples],
                                              [x[1] for x in self.outcome_tuples])

            mr_estimates = np.zeros((len(random_outcome), 3))
            for i in range(len(random_outcome)):
                mr_estimates[i, :] = self.do_single_term_mr_estimate(
                    (random_exposure[i], self.exposure_tuples[i][1]),
                    (random_outcome[i], self.outcome_tuples[i][1]))

            return random_exposure, random_outcome, mr_estimates

        def leave_one_out_residual_sum_of_squares(estimation_data,
                                                  weighted_outcome,
                                                  weighted_exposure):

            estimation_data = np.asarray(estimation_data)
            leave_one_out_ivw = np.zeros(shape=(len(estimation_data), 3))
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
        outcome = np.asarray(self.outcome_tuples)
        exposure = np.asarray(self.exposure_tuples)
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
        ivw_no_outliers = (np.nan, np.nan, np.nan)
        if local_p_val is not None and sum(local_p_val < significance_thresh):
            outliers = local_p_val < significance_thresh

            exposure_betas = [self.exposure_tuples[i][0] for i in range(len(self.estimation_data)) if not outliers[i]]
            outcome_betas = [self.outcome_tuples[i][0] for i in range(len(self.estimation_data)) if not outliers[i]]
            weights = [1 / self.outcome_tuples[i][1] ** 2 for i in range(len(self.estimation_data)) if not outliers[i]]

            outlier_corrected_ivw_result = WLS(exog=exposure_betas, endog=outcome_betas, weights=weights).fit()




        return outlier_corrected_ivw_result.params[0], outlier_corrected_ivw_result.bse[0], outlier_corrected_ivw_result.pvalues[0]



    def do_lda_mr_egger(self, ld_matrix):

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


        :param list_of_outcome_tuples:
        :param list_of_exposure_tuples:
        :param pearson_ld_matrix:
        :return:
        """


        if len(list_of_outcome_tuples) < 3:
            raise ValueError("Could not do lda simulate_mr egger on estimates, too little estimates supplied")

        marginal_exposure = np.asarray(list_of_exposure_tuples)
        marginal_outcome = np.asarray(list_of_outcome_tuples)

        #flip to make exposure strictly positive.
        to_flip = marginal_exposure[:, 0] < 0

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


##just synonym classes.

class MRPresso(MendelianRandomization):

    def __init__(self):
        super().__init__()


class LDAMREgger(MendelianRandomization):

    def __init__(self):
        super().__init__()