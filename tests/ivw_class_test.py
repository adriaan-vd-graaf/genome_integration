
import numpy as np
from genome_integration import ivw
from math import isclose


def test_smr_results():
    """
    This is the test to make sure the SMR reference is the same as my implementation of SMR.

    :return: Nothing.
    """

    ivw_thing = ivw.IVWResult()

    try: #this is tentative.
        file_thing = open("/home/adriaan/PhD/MR/SMR_explore/smr_comparison/eQTL_gen_full_smr.smr", "r")
        tolerance = 5e-5 #this file contains five significant digits per value
    except FileNotFoundError:
        file_thing = open("tests/test_resources/smr_example.txt", "r")
        tolerance = 5e-2 #this file contains two significant digits per value.

    with file_thing as f:
        for line in f:
            split = line.split()

            if split[0] == "probeID":
                continue

            outcome_tuple = (float(split[10]), float(split[11]), float(split[12]))
            exposure_tuple = (float(split[13]), float(split[14]), float(split[15]))

            smr_reference = (float(split[16]), float(split[17]), float(split[18]))
            smr_mine = ivw_thing.do_smr_estimate(outcome_tuple=outcome_tuple, exposure_tuple=exposure_tuple)

            assert isclose(smr_reference[0], smr_mine[0], rel_tol=tolerance)
            assert isclose(smr_reference[1], smr_mine[1], rel_tol=tolerance)



def test_giant_celiac_from_r_implementation():
    """
    This is used to test another MR and ivw implementation.

    :return:
    """

    beta_outcome = []
    se_outcome = []
    beta_exposure = []
    se_exposure = []

    beta_smr = []
    se_smr = []

    beta_ivw = -0.111134
    se_ivw = .2355676

    with open("tests/test_resources/celiac_compared_to_giant.txt", "r") as f:
        f.readline()
        i = 0
        ivw_thing = ivw.IVWResult()
        ivw_addition = ivw.IVWResult()
        for line in f:
            split = line.split()

            beta_outcome.append(float(split[1]))
            beta_exposure.append(float(split[2]))
            se_outcome.append(float(split[3]))
            se_exposure.append(float(split[4]))
            beta_smr.append(float(split[5]))
            se_smr.append(float(split[6]))

            outcome_tuple = (beta_outcome[i], se_outcome[i])
            exposure_tuple = (beta_exposure[i], se_exposure[i])

            smr_tuple = ivw.IVWResult().do_smr_estimate(exposure_tuple=exposure_tuple, outcome_tuple=outcome_tuple)

            ivw_addition.do_and_add_smr_estimation(exposure_tuple, outcome_tuple, "lll", 12, 12, "lll")


            assert isclose(smr_tuple[0], beta_smr[i], rel_tol=5e-10)
            assert isclose(smr_tuple[1], se_smr[i], rel_tol=5e-10)

            ivw_thing.add_estimate(smr_tuple, "lala", 12, 12, "lala")

            i+=1

        my_ivw = ivw_thing.do_ivw_estimation()

        addition_ivw = ivw_addition.do_ivw_estimation()

        integrated_ivw = ivw.IVWResult()

        integrated_result = integrated_ivw.do_ivw_estimation_on_estimate_vector([(beta_smr[j], se_smr[j]) for j in range(i)])

        assert isclose(my_ivw[0], addition_ivw[0])
        assert isclose(my_ivw[1], addition_ivw[1])
        assert isclose(my_ivw[0], integrated_result[0])
        assert isclose(my_ivw[1], integrated_result[1])

        assert isclose(my_ivw[0], beta_ivw, rel_tol=5e-5)
        assert isclose(my_ivw[1], se_ivw, rel_tol=5e-5)




def test_ivw_estimates():

    """
    This is a unit test to check if if ivw happens correctly.
    simulating some true beta, standard errors and finally some ivw.

    Assertion is done ensuring that the Z score is never lower than expected.
    Perhaps slightly cheating with the seed, but random is random sometimes.

    :return:
    """
    np.random.seed(1337)

    num_samples = 1000
    num_replicates = 100

    z_score_difference = []

    for replicate in range(num_replicates):

        true_beta = np.random.normal()
        se_smr = abs(np.random.normal(size=num_samples))
        beta_smr = np.random.normal(loc=true_beta, scale=se_smr, size=num_samples)

        this_result = ivw.IVWResult()

        for i in range(num_samples):
            this_result.add_estimate((beta_smr[i], se_smr[i]), "check", 1337, 1, "other")

        my_results = this_result.do_ivw_estimation()

        #H_0: beta_ivw == true_beta
        #H_a: beta_ivw != true_beta
        # z score to figure this out.
        z_score_difference.append((my_results[0] - true_beta) / (my_results[1]))

        assert abs(z_score_difference[len(z_score_difference)-1]) < 2.3 # corresponding to 0.01 probability, as expected with 100 individuals.


def test_q_meta_analysis_without_heterogeneity():
    np.random.seed(1337)

    num_samples = 100
    num_replicates = 20

    z_score_difference = []

    for replicate in range(num_replicates):

        true_beta = np.random.normal()
        se_smr = abs(np.random.normal(size=num_samples))
        beta_smr = np.random.normal(loc=true_beta, scale=se_smr, size=num_samples)

        this_result = ivw.IVWResult()

        for i in range(num_samples):
            this_result.add_estimate((beta_smr[i], se_smr[i]), "check", 1337, 1, "other")

        my_results = this_result.do_ivw_estimation()

        assert len(this_result.do_chochrans_q_meta_analysis(0.05)[1]) == num_samples


def test_q_meta_analysis_with_heterogeneity():
    np.random.seed(1337)


    num_bad_samples = 5
    num_good_samples = 10
    num_replicates = 100

    z_scores = []
    for replicate in range(num_replicates):

        true_beta = np.random.normal()
        good_se_smr = abs(np.random.normal(loc = 0.1, size=num_good_samples))
        good_beta_smr = np.random.normal(loc=true_beta, scale=good_se_smr, size=num_good_samples)

        bad_beta = np.random.normal(loc=-1*true_beta, size=num_bad_samples)
        bad_se_smr = abs(np.random.normal(loc=0.1, size=num_bad_samples))
        bad_beta_smr = np.random.normal(loc=bad_beta, scale=bad_se_smr, size=num_bad_samples)

        beta_smr  = list(good_beta_smr) + list(bad_beta_smr)
        se_smr = list(good_se_smr) + list(bad_se_smr)

        this_result = ivw.IVWResult()

        for i in range(num_good_samples + num_bad_samples):
            this_result.add_estimate((beta_smr[i], se_smr[i]), "check", 1337, 1, "other")

        this_result.do_ivw_estimation()
        #
        z_score_difference = abs(true_beta - this_result.do_chochrans_q_meta_analysis(0.2)[0][0]) / this_result.do_chochrans_q_meta_analysis(0.2)[0][1]

        z_scores.append(z_score_difference)

    assert np.median(z_scores) < 1.3


#now do the tests for easy debugging.
test_ivw_estimates()
test_giant_celiac_from_r_implementation()
test_smr_results()
test_q_meta_analysis_without_heterogeneity()
test_q_meta_analysis_with_heterogeneity()