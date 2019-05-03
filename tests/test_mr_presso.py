from genome_integration import causal_inference

def test_mr_presso_with_reference_implementation():
    """
    This test tests the reference implementation of MR presso, found on github:
    https://github.com/rondolab/MR-PRESSO

    It loads in the dataset located in the data folder, and then does MR presso on the two exposures from this dataset.
    on the same outcome in this dataset.

    The functions used for this are the following (adapted from reference implementation):
    mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = SummaryStats, NbDistribution = 1000,  SignifThreshold = 0.2)
    mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E2_effect", SdOutcome = "Y_se", SdExposure = "E2_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = SummaryStats, NbDistribution = 1000,  SignifThreshold = 0.2)
    :return:
    """


    #these results are from the MR presso function, which was used in
    # the results are fluid, as it is a permutation scheme, but should be very close to what I find.
    mr_presso_reference_ex1 = (0.5014829, 0.01047948)
    mr_presso_reference_ex2 = (0.8613234, 0.02143151)

    mr_presso_object_ex1 = causal_inference.MendelianRandomization()
    mr_presso_object_ex2 = causal_inference.MendelianRandomization()

    resource_path = '/'.join(('test_resources', 'mr_presso_data.txt'))

    if len(__file__.split("/")) >1:
        mr_presso_file =  "{}/{}".format("/".join(__file__.split("/")[:-1]), resource_path)
    else:
        mr_presso_file = resource_path

    with open(mr_presso_file, "r") as f:
        f.readline()
        for line in f:
            split = line.split()
            mr_presso_object_ex1.do_and_add_single_term_mr_estimation(
                (float(split[0]), float(split[1])) , (float(split[6]), float(split[7]))
            )
            mr_presso_object_ex2.do_and_add_single_term_mr_estimation(
                (float(split[3]), float(split[4])), (float(split[6]), float(split[7]))
            )

    mr_presso_result_ex1 = mr_presso_object_ex1.mr_presso(n_sims=1000, significance_thresh=0.2)
    mr_presso_result_ex2 = mr_presso_object_ex2.mr_presso(n_sims=1000, significance_thresh=0.2)


    # The precision of these results is dependent on the n_sims
    # if you want you can reduce the precision (0.02 below) and subsequently increase the n_sims parameter,
    # but this takes a fair bit of time.

    assert(mr_presso_result_ex1[0] - mr_presso_reference_ex1[0] < 0.02)
    assert (mr_presso_result_ex1[1] - mr_presso_reference_ex1[1] < 0.02)


    assert abs(mr_presso_result_ex2[0] - mr_presso_reference_ex2[0]) < 0.02
    assert abs(mr_presso_result_ex2[1] - mr_presso_reference_ex2[1]) < 0.02