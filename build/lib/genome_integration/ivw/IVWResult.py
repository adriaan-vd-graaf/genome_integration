"""
This class is used to do a 'simple'
inverse variance estimate on independent effects.
"""

__author__      = "Adriaan van der Graaf"
__copyright__   = "Copyright 2017, Adriaan van der Graaf"

import math
import scipy
import numpy as np
from .. import file_utils


class IVWResult:

    def __init__(self):
        self.estimation_done = False  # when the estimation is called, this will return true
        self.estimation_data = []     # will be filled with tuples containing betas and standard error
        self.estimation_snps = []


        #these are optional

        self.outcome_tuples = []
        self.exposure_tuples = []

        self.beta_ivw = np.nan
        self.se_ivw = np.nan
        self.p_value = np.nan

    # this is perhaps slow as you are appending betas all the time.
    # but I expect there to be at most 100 betas, so later, if it seems slow, adding the functionality.

    def add_estimate(self, beta_se_tuple, snp_name_1, pos, chr, snp_name_2):
        self.estimation_done = False
        self.estimation_data.append(beta_se_tuple)
        self.estimation_snps.append((snp_name_1, pos, chr, snp_name_2))

    def do_ivw_estimation(self):

        if len(self.estimation_data) == 0:
            raise RuntimeError('No estimates supplied to do estimation')

        beta_top = 0.0
        beta_bottom = 0.0

        for smr_result in self.estimation_data:
            beta_top += smr_result[0] * (smr_result[1] ** -2)
            beta_bottom += (smr_result[1] ** -2)

        self.beta_ivw = beta_top / beta_bottom
        self.se_ivw = math.sqrt(1 / beta_bottom)

        self.p_value = scipy.stats.norm.sf(abs(self.beta_ivw / self.se_ivw)) * 2

        self.estimation_done = True

        return self.beta_ivw, self.se_ivw, self.p_value

    def do_ivw_estimation_on_estimate_vector(self, estimation_vec):

        if len(estimation_vec) == 0:
            raise RuntimeError('No estimates supplied to do estimation')

        beta_top = 0.0
        beta_bottom = 0.0

        for smr_result in estimation_vec:
            beta_top += smr_result[0] * (smr_result[1] ** -2)
            beta_bottom += (smr_result[1] ** -2)

        beta_ivw = beta_top / beta_bottom
        se_ivw = math.sqrt(1 / beta_bottom)

        p_value = scipy.stats.norm.sf(abs(beta_ivw / se_ivw)) * 2

        return beta_ivw, se_ivw, p_value


    def get_ivw_estimates(self):
        if self.estimation_done:
            return self.beta_ivw, self.se_ivw, self.p_value
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
            self.beta_ivw,
            self.se_ivw,
            self.p_value,
            ','.join([x[0] for x in self.estimation_snps]),
            ','.join([str(x[0]) for x in self.estimation_data]),
            ','.join([str(x[1]) for x in self.estimation_data])
        )

    def write_ivw_estimates(self, filename):
        strings = ["beta_ivw\tse_ivw\tp_val_ivw\testimate_type\tsnp_1\tbp\tchr\tsnp_2",
                   "\t".join(
                            [
                            str(self.beta_ivw),
                            str(self.se_ivw),
                            str(self.p_value),
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
