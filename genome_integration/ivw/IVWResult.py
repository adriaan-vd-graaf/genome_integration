"""
This class is used to do a 'simple'
inverse variance estimate on independent effects.
"""

__author__      = "Adriaan van der Graaf"
__copyright__   = "Copyright 2017, Adriaan van der Graaf"

import math
import copy
import numpy as np
import scipy.stats
from statsmodels.regression.linear_model import WLS
from statsmodels.tools.tools import add_constant


from .. import file_utils
from .. import association


class IVWResult(association.BaseAssociation):

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



    # this is perhaps slow as you are appending betas all the time.
    # but I expect there to be at most 100 betas, so later, if it seems slow, adding the functionality.
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

        self.beta, self.se, self.wald_p_val = self.do_ivw_estimation_on_estimate_vector(self.estimation_data)
        self.estimation_done = True
        return self.beta, self.se, self.wald_p_val

    def do_ivw_estimation_on_estimate_vector(self, estimation_vec):

        if len(estimation_vec) == 0:
            raise RuntimeError('No estimates supplied to do estimation')

        beta_top = 0.0
        beta_bottom = 0.0

        i = 0
        for smr_result in estimation_vec:

            # make sure the standard error cannot be zero.
            if smr_result[1] == 0:
                smr_result[1] = np.nextafter(0,1)

            self.ivw_intermediate_top.append(smr_result[0] * (smr_result[1] ** -2))
            self.ivw_intermediate_bottom.append((smr_result[1] ** -2))

            beta_top += self.ivw_intermediate_top[i]
            beta_bottom += self.ivw_intermediate_bottom[i]
            i += 1

        beta_ivw = beta_top / beta_bottom
        se_ivw = math.sqrt(1 / beta_bottom)

        p_value = scipy.stats.norm.sf(abs(beta_ivw / se_ivw)) * 2

        return beta_ivw, se_ivw, p_value


    def get_ivw_estimates(self):
        if self.estimation_done:
            return self.beta, self.se, self.wald_p_val
        else:
            raise RuntimeError('The estimation is not done, so no estimates were made.')

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

        return beta_smr, se_smr, p_value

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

        file_utils.write_list_to_newline_separated_file(strings, filename)


    def do_chochrans_q_meta_analysis(self, p_value_threshold):

        if len(self.estimation_data) < 3:
            raise ValueError("Less than three estimates supplied, cannot do cochrans q analysis")

        if not self.estimation_done:
            raise ValueError("No ivw estimation done.")


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

            #determine the new ivw estimates
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

        num_estimates = len(self.estimation_data)

        # runtime checks.

        if num_estimates < 4:
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