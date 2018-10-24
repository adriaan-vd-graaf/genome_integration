from .. import association
from .. import ivw
import numpy as np
import scipy.stats


class MRPresso(ivw.IVWResult):

    def __init__(self):
        super().__init__()


    def do_single_term_mr_estimate(self, exposure_tuple, outcome_tuple):
        """
        Determine SMR test effect and standard error.

        :param exposure_data: Exposure estimates which must have the methods get_beta and get_z_score
        :param outcome_data: Outcome estimates which must have the methods get_beta and get_z_score
        :return: tuple of smr beta and smr se
        """

        beta_mr = outcome_tuple[0] / exposure_tuple[0]
        se_mr = np.sqrt((outcome_tuple[1] **2)  / (exposure_tuple[0] **2))

        z_score = beta_mr / se_mr
        p_value = scipy.stats.norm.sf(abs(z_score)) * 2

        return [beta_mr, se_mr, p_value]


    def do_and_add_single_term_mr_estimation(self, exposure_tuple, outcome_tuple):
        estimates = self.do_single_term_mr_estimate(exposure_tuple, outcome_tuple)
        self.estimation_data.append(estimates)

        self.outcome_tuples.append(outcome_tuple)
        self.exposure_tuples.append(exposure_tuple)



    def mr_presso(self, n_sims = 1000, significance_thresh=0.05):


        def make_random_data():
            beta_ivw, _, _ =  self.get_ivw_estimates()
            random_exposure = np.random.normal([x[0] for x in self.exposure_tuples], [x[1] for x in self.exposure_tuples])
            random_outcome = np.random.normal([beta_ivw * x[0] for x in self.exposure_tuples],
                                              [x[1] for x in self.outcome_tuples])

            mr_estimates = np.zeros((len(random_outcome), 3))
            for i in range(len(random_outcome)):
                mr_estimates[i,:] = self.do_single_term_mr_estimate( (random_exposure[i], self.exposure_tuples[i][1]),
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

            rss = (weighted_outcome - leave_one_out_ivw[:,0] * weighted_exposure) ** 2

            return rss, leave_one_out_ivw

        def make_random_data_and_return_rss(weights):
            exposure, outcome, mr_estimates = make_random_data()

            weighted_exposure = exposure * weights
            weighted_outcome = outcome * weights

            rss, _ = leave_one_out_residual_sum_of_squares(mr_estimates, weighted_outcome, weighted_exposure)

            return np.sum(rss), np.concatenate((exposure.reshape(len(exposure), 1), outcome.reshape(len(outcome), 1)), axis=1)

        def randomly_sample_distortion(outlier_indices):
            estimates = np.asarray(self.estimation_data)
            estimates_no_outliers = np.delete(estimates, outlier_indices, axis=0)
            estimates_only_outliers = estimates[outlier_indices,:][0]

            indices_sampled_from_no_outliers = np.random.choice(estimates_no_outliers.shape[0],
                                                                size=estimates_no_outliers.shape[0],
                                                                replace=True)
            return self.do_ivw_estimation_on_estimate_vector(
                np.concatenate((estimates_no_outliers[indices_sampled_from_no_outliers,:], estimates_only_outliers))
            )

        # runtime checks.
        num_estimates = len(self.estimation_data)

        if num_estimates < 4:
            raise ValueError("Only {} estimates supplied, need at least three to find simulate_mr presso outliers".format(num_estimates))

        if len(self.exposure_tuples) != num_estimates:
            raise ValueError("No exposure sumstats present, cannot do mr_presso outlier.")

        if len(self.outcome_tuples) != num_estimates:
            raise ValueError("No outcome sumstats present, cannot do mr_presso outlier.")


        #this is just following MR presso.
        outcome = np.asarray(self.outcome_tuples)
        exposure = np.asarray(self.exposure_tuples)
        weighted_outcome = np.asarray([x[0] / np.sqrt(x[1]**2) for x in self.outcome_tuples], dtype=float)
        weighted_exposure = np.asarray([self.exposure_tuples[i][0] / np.sqrt(self.outcome_tuples[i][1]**2)
                                        for i in range(len(self.exposure_tuples))], dtype=float)
        weights = np.asarray([1 / np.sqrt(x[1] **2) for x in self.outcome_tuples], dtype=float)

        rss, list_of_assocs = leave_one_out_residual_sum_of_squares(self.estimation_data,
                                                                    weighted_outcome,
                                                                    weighted_exposure)


        expected_results = [make_random_data_and_return_rss(weights) for _ in range(n_sims)]

        sim_rss = [x[0] for x in expected_results]

        global_p_val = sum( sim_rss > sum(rss) ) / n_sims
        local_p_val = None
        if global_p_val < significance_thresh:
            expected_betas = np.zeros((num_estimates, n_sims, 2), dtype=float)
            for i in range(n_sims):
                expected_betas[:,i] = expected_results[i][1]

            difference = outcome[:,0] - exposure[:,0] * list_of_assocs[:,0]
            expected_difference = expected_betas[:,:,1] - expected_betas[:,:,0] * np.tile(list_of_assocs[:,0], (n_sims, 1)).transpose()
            local_p_val = np.sum(expected_difference ** 2 > (difference ** 2).reshape((len(difference),1)), axis = 1) / n_sims
            local_p_val = np.asarray([x * len(difference) if x * len(difference) < 1.0 else 1.0  for x in local_p_val])

        #distortion test.
        ivw_no_outliers = (np.nan, np.nan, np.nan)
        if local_p_val is not None and sum(local_p_val < significance_thresh):
            outliers = np.where(local_p_val < significance_thresh)
            ivw_all = self.get_ivw_estimates()
            ivw_no_outliers  = self.do_ivw_estimation_on_estimate_vector(
                np.delete(np.asarray(self.estimation_data), outliers, axis=0)
            )
            observed_bias = ivw_all[0]

            expected_ivw = np.asarray([randomly_sample_distortion(outliers) for _ in range(n_sims)])


        return ivw_no_outliers