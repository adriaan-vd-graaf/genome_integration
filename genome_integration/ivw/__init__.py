"""
This class is used to do a 'simple'
inverse variance estimate on independent effects.
"""

__author__      = "Adriaan van der Graaf"
__copyright__   = "Copyright 2017, Adriaan van der Graaf"

import math
import scipy
from .. import file_utils


class IVWResult:

    def __init__(self):
        self.estimation_done = False  # when the estimation is called, this will return true
        self.estimation_data = []     # will be filled with tuples containing betas and standard error
        self.estimation_snps = []

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

        self.00 = beta_top / beta_bottom
        self.se_ivw = math.sqrt(1 / beta_bottom)

        self.p_value = scipy.stats.norm.sf(abs(self.beta_ivw / self.se_ivw)) * 2

        self.estimation_done = True

        return self.beta_ivw, self.se_ivw, self.p_value

    def get_ivw_estimates(self):
        if self.estimation_done:
            return self.beta_ivw, self.se_ivw, self.p_value
        else:
            raise RuntimeError('The estimation is not done, so no estimates were made.')

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
