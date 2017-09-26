
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
        file_thing = open("test_resources/smr_example.txt", "r")
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



    # weirdly this 'seems' to be biassed towards lower values, perhaps meaning that we're underestimating the effect.
    # print(np.mean(z_score_difference))
    # print(np.sum(np.sign(z_score_difference)))

    # print(np.max(np.abs(z_score_difference)))




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


#this does not work... so we're going to simulate beta_smr, and then make
def test_smr_ivw_integration():
    """
    This script will make sure that there

    """
    np.random.seed(1337)

    num_samples = 100
    num_replicates = 1

    z_score_difference = []

    for replicate in range(num_replicates):

        true_beta = np.random.normal()
        se_smr = abs(np.random.normal(size=num_samples))
        beta_smr = np.random.normal(loc=true_beta, scale=se_smr, size=num_samples)

        z_sq_smr = beta_smr ** 2 / se_smr ** 2

        """
        
        Now we do some math. We know the following:
        
        beta_smr^2  / se_smr^2 = (z_outcome^2 * z_exposure^2) / (z_outcome^2 + z_exposure^2)
        
        0
        after some algebra, I come to the following:
        
        z_outcome^2 = - ( (beta_smr^2  / se_smr^2 ) * z_exposure^2 ) / ( (beta_smr^2  / se_smr^2) - z_exposure^2 ))
         
        se_outcome = (beta_outcome * np.sqrt( (beta_exposure^2 / (beta_smr^2  / se_smr^2 ) - se_exposure^2)) / beta_exposure 
        
        If I then simulate some beta_exposure, with a se_exposure, 
        I have z_exposure, and through some algebra, I get z_outcome as well.
        
        """

        se_exposure = np.abs(np.random.normal(size=num_samples))
        beta_exposure = np.random.normal(scale=se_exposure, size=num_samples)

        z_exposure = beta_exposure / se_exposure
        z_outcome = np.sqrt(np.abs(( z_sq_smr * z_exposure ** 2) / ( z_sq_smr - z_exposure ** 2)))


        se_outcome = np.abs(( beta_smr * beta_exposure  * np.sqrt(np.abs(beta_exposure ** 2 / (beta_smr ** 2 / se_smr ** 2) - se_exposure ** 2))) / beta_exposure)

        print(se_outcome)

        beta_outcome = np.random.normal(z_outcome * se_outcome, se_outcome, num_samples)

        print(np.mean(beta_outcome / beta_exposure))

        print(true_beta)


test_smr_ivw_integration()
test_ivw_estimates()
test_smr_results()
