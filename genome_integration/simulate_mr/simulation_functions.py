import copy
import numpy as np
import statsmodels.api as sm
from . import geno_functions


def simulate_phenotypes_also_return_unobserved_phenotype(exposure_1_causal, exposure_2_causal,
                                                            exposure_1_n_causal, exposure_2_n_causal,
                                                            exposure_geno, exposure_ld,
                                                            outcome_geno,
                                                            overlapping_causal_snps=0,
                                                            error_sd=1.0,
                                                            confounder_sd=0.5,
                                                            inside_phi=0.0,
                                                            directional_pleiotropy=True
                                                            ):
    """
    This function simulates two cohorts which are under the same genetic control
    This function is present to ensure that across simulation scenarios, the same parameters are used.

    Step by step:

    1.
        exposure 1 SNPs are chosen, and they are chosen not to be in higher ld than 0.95 R ^ 2

    2.
        exposure 2 SNPs are chosen in two ways. `overlapping_causal_snps` are chosen without replacement from the
        exposure 1 SNPs, the rest are chosen from SNPs that are in LD with exposure 1.

    3.
        The confounder is simulated, based on some error and if `inside_phi` is nonzero, it is given a genetic effect.

    4.
        The effect sizes of exposure 1 and of exposure 2 eQTLs are chosen.

    5.
        The phenotypes are simulated based on the function simulate_phenotypes_extended,
        in the exposure cohort, and the outcome cohort indepdently.

    6.
        returning the following:
            exposure_phenotype, outcome_phenotype,
            exposure_1_causal_snps, exposure_1_betas,
            exposure_2_causal_snps, exposure_2_betas



    :param exposure_1_causal:
    The causal effect on exposure 1

    :param exposure_2_causal:
    The causal effect on exposure 2

    :param exposure_1_n_causal:
    The number of causal SNPs for exposure 1e

    :param exposure_2_n_causal:
    The number of causal SNPs for exposure 2

    :param exposure_geno:
    The genotype matrix for the exposure

    :param exposure_ld:
    The ld R^2 matrix for the exposure genotypes

    :param outcome_geno:
    The genotype matrix for the outcome

    :param outcome_ld:
    The ld R^2 matrix for the outcome genotypes. currently not used.

    :param overlapping_causal_snps:
    An integer number representing the number of causal SNPs for exposure 2 that overlap exposure 1.

    :param error_sd:
    The scale (standard deviation) parameter of the error that is added to a phenotype

    :param inside phi:
    The phi parameter (from the egger paper notation), that makes the confounder dependent of the genotype.
    if non zero, a genetic effect with effect size (0,2) is added to the confounder, using the exposure 2 causal SNPs

    :param directional_pleiotropy:
    if there should be directional pleiotropy or not.

    :return:
    Returns a tuple of the following numpy arrays:

    1 by n vector of simulated phenotypes for the exposure,
    1 by n vector simulated phenotypes for the outcome,
    1 by n causal vector of the indices of the causal snps of exposure 1.
    1 by n causal vector of the effect sizes of the causal snps of exposure 1.
    1 by n causal vector of the indices of the causal snps of exposure 2.
    1 by n causal vector of the effect sizes of the causal snps of exposure 2.

    """

    # runtime check
    if overlapping_causal_snps > exposure_1_n_causal:
        raise RuntimeError("The total number of overlapping snps cannot be larger than the number of causal snps for " +
                           "exposure 1 ")
    if overlapping_causal_snps > exposure_2_n_causal:
        raise RuntimeError("The total number of overlapping snps cannot be larger than the number of causal snps for " +
                           "exposure 2 ")

    exposure_beta_limits = (-.5, .5)  # the limits of the normal distribution for the exposure.

    upper_ld_bound = 0.95  # The upper r^2 boundary for causal SNP selection of linkage
    lower_ld_bound = 0.25  # The lower r^2 boundary for cauasal SNP selection of linkage

    exposure_ld = exposure_ld ** 2

    """
    Step 1.
    """

    # make sure the exposure Causal SNPs are not in > upper_ld_bound R^2 with another.
    exposure_1_causal_snps = np.zeros([exposure_1_n_causal], dtype=int)
    i = 0
    exposure_1_causal_snps[i] = np.random.choice(exposure_geno.shape[0], 1, replace=False)
    i += 1

    while i < exposure_1_n_causal:
        # choosing SNPS that are not in full linkage (R**2 > 0.95) with one another, but do share some linkage.
        lower_than_upper = np.sum(exposure_ld[exposure_1_causal_snps[:i], :] < upper_ld_bound, axis=0) == i
        choice = np.random.choice(np.where(lower_than_upper)[0], 1, replace=False)
        exposure_1_causal_snps[i] = choice
        i += 1

    """
    Step 2.
        Select overlapping.
        If there are non overlapping SNPs left, select SNPs in LD with exposure 1.
    """
    # Here, the snps of exposure are overlapped, depending on the `overlapping_causal_snps` number, could be zero.
    exposure_2_overlapping_snps = np.random.choice(exposure_1_causal_snps, overlapping_causal_snps, replace=False)

    """
        Here we select SNPs for exposure 2, that are in LD with exposure 1.
        But we do not look for LD with already overlapping causal snps for exposure 1.

        Only if there are no overlapping snps left for exposure 2 will we just sat our exposure 2 causal snps are done.
    """
    if exposure_2_n_causal - overlapping_causal_snps > 0:

        exposure_1_causal_snps_non_overlapping = np.asarray(
            [x for x in exposure_1_causal_snps if x not in exposure_2_overlapping_snps], dtype=int)

        # permute this vector, so ordering is random.
        permuted_exposure_1_causal_non_overlapping = np.random.permutation(exposure_1_causal_snps_non_overlapping)
        exposure_2_ld_snps = np.zeros((exposure_2_n_causal - overlapping_causal_snps), dtype=int)

        # iterate over all causal snps, and choose one that is within the ld window
        i = 0
        while i < (exposure_2_n_causal - overlapping_causal_snps):

            ld_with_causal = exposure_ld[permuted_exposure_1_causal_non_overlapping[i], :]

            try:
                chosen_indice = np.random.choice(
                    np.where((ld_with_causal > lower_ld_bound) & (ld_with_causal < upper_ld_bound))[0], 1)
                exposure_2_ld_snps[i] = chosen_indice
                i += 1
            except:
                raise ValueError("Could not find SNPs in an ld window around the exposure snp.")

        exposure_2_causal_snps = np.concatenate((exposure_2_overlapping_snps, exposure_2_ld_snps))

    else:
        exposure_2_causal_snps = exposure_2_overlapping_snps

    """
    Step 3
        Confounder, including the inside assumption
    """

    confounder_exposure_cohort = np.random.normal(0, confounder_sd, exposure_geno.shape[1])
    confounder_outcome_cohort = np.random.normal(0, confounder_sd, outcome_geno.shape[1])

    if inside_phi != 0.0:
        phi_values = np.random.uniform(0, inside_phi, (1, len(exposure_2_causal_snps)))
        confounder_exposure_cohort += (phi_values @ exposure_geno[exposure_2_causal_snps, :]).reshape(
            exposure_geno.shape[1])
        confounder_outcome_cohort += (phi_values @ outcome_geno[exposure_2_causal_snps, :]).reshape(
            outcome_geno.shape[1])

    """
    Step 4
        Simulate the effect sizes of the causal snps
    """

    exposure_1_betas = np.random.uniform(exposure_beta_limits[0], exposure_beta_limits[1],
                                         exposure_1_n_causal)
    if directional_pleiotropy:
        exposure_2_betas = np.random.uniform(0.0, exposure_beta_limits[1],
                                             exposure_2_n_causal)
    else:
        exposure_2_betas = np.random.uniform(exposure_beta_limits[0], exposure_beta_limits[1],
                                             exposure_2_n_causal)

    """
    Step 5
        The actual phenotypes are simulated.
    """

    # exposure 1 is of interest. after simulation of phenotypes, exposure 2 is discarded.
    exposure_phenotype, exposure_2_phenotype, _ = simulate_phenotypes_extended(exposure_geno,
                                                            exposure_1_causal_snps=exposure_1_causal_snps,
                                                            exposure_1_causal_effect=exposure_1_causal,
                                                            exposure_1_betas=exposure_1_betas,
                                                            exposure_2_causal_snps=exposure_2_causal_snps,
                                                            exposure_2_causal_effect=exposure_2_causal,
                                                            exposure_2_betas=exposure_2_betas,
                                                            error_scale=error_sd,
                                                            confounder=confounder_exposure_cohort
                                                            )
    # center to zero, set variance to 1
    exposure_phenotype -= np.mean(exposure_phenotype)
    exposure_phenotype /= np.std(exposure_phenotype)

    _, _, outcome_phenotype = simulate_phenotypes_extended(outcome_geno,
                                                           exposure_1_causal_snps=exposure_1_causal_snps,
                                                           exposure_1_causal_effect=exposure_1_causal,
                                                           exposure_1_betas=exposure_1_betas,
                                                           exposure_2_causal_snps=exposure_2_causal_snps,
                                                           exposure_2_causal_effect=exposure_2_causal,
                                                           exposure_2_betas=exposure_2_betas,
                                                           error_scale=error_sd,
                                                           confounder=confounder_outcome_cohort
                                                           )

    # center to zero, set variance to 1
    outcome_phenotype -= np.mean(outcome_phenotype)
    outcome_phenotype /= np.std(outcome_phenotype)

    return exposure_phenotype, outcome_phenotype, \
           exposure_1_causal_snps, exposure_1_betas, \
           exposure_2_causal_snps, exposure_2_betas, \
           exposure_2_phenotype




def simulate_phenotypes(exposure_1_causal, exposure_2_causal,
                        exposure_1_n_causal, exposure_2_n_causal,
                        exposure_geno, exposure_ld,
                        outcome_geno,
                        overlapping_causal_snps = 0,
                        error_sd = 1.0,
                        confounder_sd = 0.5,
                        inside_phi = 0.0,
                        directional_pleiotropy = True
                        ):
    """
    This function simulates two cohorts which are under the same genetic control
    This function is present to ensure that across simulation scenarios, the same parameters are used.

    Step by step:

    1.
        exposure 1 SNPs are chosen, and they are chosen not to be in higher ld than 0.95 R ^ 2

    2.
        exposure 2 SNPs are chosen in two ways. `overlapping_causal_snps` are chosen without replacement from the
        exposure 1 SNPs, the rest are chosen from SNPs that are in LD with exposure 1.

    3.
        The confounder is simulated, based on some error and if `inside_phi` is nonzero, it is given a genetic effect.

    4.
        The effect sizes of exposure 1 and of exposure 2 eQTLs are chosen.

    5.
        The phenotypes are simulated based on the function simulate_phenotypes_extended,
        in the exposure cohort, and the outcome cohort indepdently.

    6.
        returning the following:
            exposure_phenotype, outcome_phenotype,
            exposure_1_causal_snps, exposure_1_betas,
            exposure_2_causal_snps, exposure_2_betas



    :param exposure_1_causal:
    The causal effect on exposure 1

    :param exposure_2_causal:
    The causal effect on exposure 2

    :param exposure_1_n_causal:
    The number of causal SNPs for exposure 1e

    :param exposure_2_n_causal:
    The number of causal SNPs for exposure 2

    :param exposure_geno:
    The genotype matrix for the exposure

    :param exposure_ld:
    The ld R^2 matrix for the exposure genotypes

    :param outcome_geno:
    The genotype matrix for the outcome

    :param outcome_ld:
    The ld R^2 matrix for the outcome genotypes. currently not used.

    :param overlapping_causal_snps:
    An integer number representing the number of causal SNPs for exposure 2 that overlap exposure 1.

    :param error_sd:
    The scale (standard deviation) parameter of the error that is added to a phenotype

    :param inside phi:
    The phi parameter (from the egger paper notation), that makes the confounder dependent of the genotype.
    if non zero, a genetic effect with effect size (0,2) is added to the confounder, using the exposure 2 causal SNPs

    :param directional_pleiotropy:
    if there should be directional pleiotropy or not.

    :return:
    Returns a tuple of the following numpy arrays:

    1 by n vector of simulated phenotypes for the exposure,
    1 by n vector simulated phenotypes for the outcome,
    1 by n causal vector of the indices of the causal snps of exposure 1.
    1 by n causal vector of the effect sizes of the causal snps of exposure 1.
    1 by n causal vector of the indices of the causal snps of exposure 2.
    1 by n causal vector of the effect sizes of the causal snps of exposure 2.

    """

    #runtime check
    if overlapping_causal_snps > exposure_1_n_causal:
        raise RuntimeError("The total number of overlapping snps cannot be larger than the number of causal snps for " +
                           "exposure 1 ")
    if overlapping_causal_snps > exposure_2_n_causal:
        raise RuntimeError("The total number of overlapping snps cannot be larger than the number of causal snps for " +
                           "exposure 2 ")

    exposure_beta_limits = (-.5, .5) # the limits of the normal distribution for the exposure.

    upper_ld_bound = 0.95 # The upper r^2 boundary for causal SNP selection of linkage
    lower_ld_bound = 0.25 # The lower r^2 boundary for cauasal SNP selection of linkage

    exposure_ld = exposure_ld ** 2

    """
    Step 1.
    """

    #make sure the exposure Causal SNPs are not in > upper_ld_bound R^2 with another.
    exposure_1_causal_snps = np.zeros([exposure_1_n_causal], dtype=int)
    i=0
    exposure_1_causal_snps[i] = np.random.choice(exposure_geno.shape[0], 1, replace=False)
    i += 1

    while i < exposure_1_n_causal:
        #choosing SNPS that are not in full linkage (R**2 > 0.95) with one another, but do share some linkage.
        lower_than_upper = np.sum(exposure_ld[exposure_1_causal_snps[:i], :] < upper_ld_bound, axis=0) == i
        choice = np.random.choice(np.where(lower_than_upper)[0], 1, replace=False)
        exposure_1_causal_snps[i] = choice
        i+=1

    """
    Step 2.
        Select overlapping.
        If there are non overlapping SNPs left, select SNPs in LD with exposure 1.
    """
    #Here, the snps of exposure are overlapped, depending on the `overlapping_causal_snps` number, could be zero.
    exposure_2_overlapping_snps = np.random.choice(exposure_1_causal_snps, overlapping_causal_snps, replace=False)

    """
        Here we select SNPs for exposure 2, that are in LD with exposure 1.
        But we do not look for LD with already overlapping causal snps for exposure 1.
        
        Only if there are no overlapping snps left for exposure 2 will we just sat our exposure 2 causal snps are done.
    """
    if exposure_2_n_causal - overlapping_causal_snps > 0:

        exposure_1_causal_snps_non_overlapping = np.asarray([x for x in exposure_1_causal_snps if x not in exposure_2_overlapping_snps], dtype=int)

        #permute this vector, so ordering is random.
        permuted_exposure_1_causal_non_overlapping = np.random.permutation(exposure_1_causal_snps_non_overlapping)
        exposure_2_ld_snps = np.zeros((exposure_2_n_causal - overlapping_causal_snps), dtype=int)

        #iterate over all causal snps, and choose one that is within the ld window
        i = 0
        while i < (exposure_2_n_causal - overlapping_causal_snps):

            ld_with_causal = exposure_ld[permuted_exposure_1_causal_non_overlapping[i],:]

            try:
                chosen_indice = np.random.choice(np.where((ld_with_causal > lower_ld_bound) & (ld_with_causal < upper_ld_bound))[0], 1)
                exposure_2_ld_snps[i] = chosen_indice
                i+=1
            except:
                raise ValueError("Could not find SNPs in an ld window around the exposure snp.")


        exposure_2_causal_snps = np.concatenate((exposure_2_overlapping_snps, exposure_2_ld_snps))

    else:
        exposure_2_causal_snps = exposure_2_overlapping_snps


    """
    Step 3
        Confounder, including the inside assumption
    """

    confounder_exposure_cohort = np.random.normal(0, confounder_sd, exposure_geno.shape[1])
    confounder_outcome_cohort = np.random.normal(0, confounder_sd, outcome_geno.shape[1])

    if inside_phi != 0.0:
        phi_values = np.random.uniform(0, inside_phi, (1, len(exposure_2_causal_snps)))
        confounder_exposure_cohort += (phi_values @ exposure_geno[exposure_2_causal_snps,:]).reshape(exposure_geno.shape[1])
        confounder_outcome_cohort += (phi_values @ outcome_geno[exposure_2_causal_snps, :]).reshape(outcome_geno.shape[1])


    """
    Step 4
        Simulate the effect sizes of the causal snps
    """

    exposure_1_betas = np.random.uniform(exposure_beta_limits[0], exposure_beta_limits[1],
                                         exposure_1_n_causal)
    if directional_pleiotropy:
        exposure_2_betas = np.random.uniform(0.0, exposure_beta_limits[1],
                                             exposure_2_n_causal)
    else:
        exposure_2_betas = np.random.uniform(exposure_beta_limits[0], exposure_beta_limits[1],
                                         exposure_2_n_causal)

    """
    Step 5
        The actual phenotypes are simulated.
    """


    # exposure 1 is of interest. after simulation of phenotypes, exposure 2 is discarded.
    exposure_phenotype, _, _ = simulate_phenotypes_extended(exposure_geno,
                                                   exposure_1_causal_snps=exposure_1_causal_snps,
                                                   exposure_1_causal_effect=exposure_1_causal,
                                                   exposure_1_betas=exposure_1_betas,
                                                   exposure_2_causal_snps=exposure_2_causal_snps,
                                                   exposure_2_causal_effect=exposure_2_causal,
                                                   exposure_2_betas=exposure_2_betas,
                                                   error_scale=error_sd,
                                                   confounder=confounder_exposure_cohort
                                                   )
    # center to zero, set variance to 1
    exposure_phenotype -= np.mean(exposure_phenotype)
    exposure_phenotype /= np.std(exposure_phenotype)

    _, _, outcome_phenotype = simulate_phenotypes_extended(outcome_geno,
                                                  exposure_1_causal_snps=exposure_1_causal_snps,
                                                  exposure_1_causal_effect=exposure_1_causal,
                                                  exposure_1_betas=exposure_1_betas,
                                                  exposure_2_causal_snps=exposure_2_causal_snps,
                                                  exposure_2_causal_effect=exposure_2_causal,
                                                  exposure_2_betas=exposure_2_betas,
                                                  error_scale=error_sd,
                                                  confounder=confounder_outcome_cohort
                                                  )

    # center to zero, set variance to 1
    outcome_phenotype -= np.mean(outcome_phenotype)
    outcome_phenotype /= np.std(outcome_phenotype)

    return exposure_phenotype, outcome_phenotype, \
           exposure_1_causal_snps, exposure_1_betas, \
           exposure_2_causal_snps, exposure_2_betas



def simulate_phenotypes_extended(geno,
                                exposure_1_causal_snps,
                                exposure_1_causal_effect,
                                exposure_1_betas,
                                exposure_2_causal_snps,
                                exposure_2_causal_effect,
                                exposure_2_betas,
                                error_scale,
                                confounder
                                ):
    """
    This function will simulate phenotypes based on parameters passed through it.
    Simulation is done through up two exposures, which can both be causal for some outcome.

    The math behind it is the following:
    exposure 1: y_{e1} = X_{1} \cdot b_{e1} + C + e
    exposure 2: y_{e2} = X_{2} \cdot b_{e2} + C + e

    outcome:  y_{o} =  y_{e1} b_{1} + X_{1} \cdot b_{p,1} + X_{2} \cdot b_{p,1} + y_{e2} b_{2} + C + e

    :param geno:
    A genotype matrix, and n individuals by m variants,
    numpy integer matrix,
    representing the genotypes of individuals.

    :param exposure_1_causal_snps:
    A vector of length i (<= m), number of causal snps for exposure 1, with indices in the genotype matrix from which to
    choose causal snps.
    X_{1} is the genotype matrix, with only these cauasal indices.

    :param exposure_1_causal_effect:
    A float representing the causal effect of exposure 1 on the outcome

    :param exposure_1_betas:
    b_{1} is a i by 1 vector of floats, representing the causal effect of the snp on the exposure.

    :param exposure_2_causal_snps:
    A vector of length j (<= m), number of causal snps for exposure 2, with indices in the genotype matrix from which to
    choose causal snps.
    X_{2} is the genotype matrix, with only these causal indices.

    :param exposure_2_causal_effect:
    A float representing the causal effect of exposure 1 on the outcome

    :param exposure_2_betas:
    b_{2} is a j by 1 vector of floats, representing the causal effect of the snp on the exposure.

    :param error_scale:
    e ~ N(0, error_scale)
    The scale parameter of the e term, drawn from a normal distribution with this standard deviation

    :param confounder:
    a n by 1 vector of confounder values per individual.

    :return:
    three length n vectors representing in order:
    phenotypes for exposure 1, exposure 2 and the outcome.

    """

    scaled_geno = np.apply_along_axis(geno_functions.scale_geno_vec, 1, geno)

    # exposure 1: y_{e1} = X_{1}b_{e1} + C + e
    exposure_1_phenotype = exposure_1_betas @ scaled_geno[exposure_1_causal_snps, :]

    exposure_1_phenotype += confounder
    exposure_1_phenotype += np.random.normal(0, error_scale, scaled_geno.shape[1])  # add some measurement error

    # exposure 2: y_{e2} = X_{2}b_{e2} + C + e
    exposure_2_phenotype = exposure_2_betas @ scaled_geno[exposure_2_causal_snps, :]
    exposure_2_phenotype += confounder
    exposure_2_phenotype += np.random.normal(0, error_scale, scaled_geno.shape[1])  # add some measurement error

    # y_{o} = y_{e1} b_{1} + X_{1}b_{p,1} + X_{2}_b_{p,1} + y_{e2} b_{2} + C + e
    outcome_phenotype = exposure_1_phenotype * exposure_1_causal_effect
    outcome_phenotype += exposure_2_phenotype * exposure_2_causal_effect

    outcome_phenotype += confounder
    outcome_phenotype += np.random.normal(0, error_scale, scaled_geno.shape[1])  # add some measurement error

    """
    Phenotypes are now simulated  
    """

    return exposure_1_phenotype, exposure_2_phenotype, outcome_phenotype



def select_random_snps_around_causal(ld_mat, causal_indices, num_extra=0, upper_threshold=0.5):
    """
    This function selects SNPs randomly around already identified SNPs, based on some LD threshold (upper_threshold).
    All returned indices are in LD lower than upper_threshold with one another.

    :param ld_mat:
    A J by J numpy matrix, representing some correlation metric between variants, where J is the number of variants,
    :param causal_indices:
    A N by 1 numpy matrix, where N is the number of variants selected to find the high ld snps around
    :param num_extra:
    Number of extra variants to look for
    :return:
    returns a N + num_extra by 1 numpy matrix, or less (depending on the ld matrix)

    """

    """
    Stage 1, Find SNPs in LD around the already identified causal indices. 
    """

    random_indices = []
    for causal_snp in causal_indices:

        causal_and_random = np.asarray(np.concatenate(([causal_snp],random_indices)), dtype=int)

        below_thresh = np.all(ld_mat[causal_and_random,:] <= upper_threshold, axis=0)

        if np.sum(below_thresh) == 0:
            raise RuntimeWarning("Could not find any SNPs in LD below threshold.")

        max_ld_of_candidates = np.max(ld_mat[causal_snp,below_thresh])

        max_ld_indice = np.where(ld_mat[causal_snp,:] == max_ld_of_candidates)

        random_indices.append(max_ld_indice[0][0])

    """
    Stage 2, find some extra SNPs around randomly chosen causal indices.
    """

    for extra_snp in range(num_extra):

        causal_for_reference = int(np.random.choice(causal_indices, 1, replace=True))
        causal_and_random = np.asarray(np.concatenate(([causal_for_reference], random_indices)), dtype=int)

        below_thresh = np.all(ld_mat[causal_and_random, :] <= upper_threshold, axis=0)
        if np.sum(below_thresh) == 0:
            raise RuntimeWarning("Could not find any SNPs in LD below threshold. Returning the SNPs that have been chosen")


        max_ld_of_candidates = np.max(ld_mat[causal_for_reference, below_thresh])

        max_ld_indice = np.where(ld_mat[causal_for_reference, :] == max_ld_of_candidates)

        random_indices.append(max_ld_indice[0][0])

    return random_indices


def randomly_choose_snps_and_flip_alleles_identify_causal(outcome_geno,
                                                          r_sq_mat,
                                                          exposure_betas,
                                                          causal_exposure_indices,
                                                          outcome_phenotype,
                                                          r_sq_threshold=0.3,
                                                          extra_random_snps=0,
                                                          n_iterations=50
                                                          ):

    causal_estimates = np.zeros((n_iterations, 4), dtype=float)
    geno = copy.deepcopy(outcome_geno.transpose()) #probably slow
    betas = copy.deepcopy(exposure_betas) #probably slow

    for i in range(n_iterations):

        random_indices = select_random_snps_around_causal(r_sq_mat, causal_exposure_indices,
                                                          upper_threshold=r_sq_threshold,
                                                          num_extra=extra_random_snps)

        design_mat = np.zeros( (geno.shape[0], len(random_indices) + 1 ), dtype=float)
        design_mat[:,0] = geno[:,causal_exposure_indices] @ betas
        design_mat[:,1:(len(random_indices)+1)] = geno[:,random_indices]

        try:

            model = sm.OLS(endog=outcome_phenotype, exog=design_mat)

            ols_fit = model.fit()

            causal_estimates[i, :] = ols_fit.params[0], ols_fit.bse[0], ols_fit.pvalues[0], ols_fit.bic

        except Exception as x:
            print("Did not find estimate for {}, with this exception: {}".format(i, x))


    return causal_estimates



def do_gcta_cojo_conditional(reference_geno, associations,  indices_of_interest, ref_loci):
    """
    Recreates the conditional analysis of the GCTA cojo paper.
    Does NOT do the stepwise regression.

    See the COJO paper (JIAN YANG 2011, Nature Genetics) supplemental methods for the derivation.

    :param reference_geno:
    :param marginal_estimates:
    :param marginal_maf:
    :return:
    """
    sample_size = associations[ref_loci[0].name].n_observations

    beta_se_maf = np.asarray(
                        [[associations[x.name].beta,
                             associations[x.name].se,
                             associations[x.name].minor_allele_frequency
                             ]
                         for x in ref_loci]
                        , dtype=float)

    full_d = 2*beta_se_maf[:,2] * (1 - beta_se_maf[:,2]) * sample_size
    s_squared = beta_se_maf[:,1] ** 2

    y_t_y = full_d * s_squared * (sample_size - 1) + full_d * beta_se_maf[:,0] **2
    median_y_t_y = np.median(y_t_y)

    adjusted_sample_size = (median_y_t_y / (full_d  * s_squared)) - (beta_se_maf[:,0] **2 / s_squared) + 1

    new_b = np.zeros( (len(indices_of_interest), len(indices_of_interest)), dtype=float)
    new_d = np.zeros((len(indices_of_interest), len(indices_of_interest)), dtype=float)

    w_sums = np.sum(reference_geno ** 2, axis=0)

    # now fill in the b and d vectors.
    for j in range(len(indices_of_interest)):
        new_d[j,j] = 2 * beta_se_maf[indices_of_interest[j],2] * \
                    (1 - beta_se_maf[indices_of_interest[j],2]) * \
                    adjusted_sample_size[indices_of_interest[j]]

        n_vec = np.asarray([x
                            if adjusted_sample_size[indices_of_interest[j]] > x
                            else adjusted_sample_size[indices_of_interest[j]]
                            for x in adjusted_sample_size[indices_of_interest]], dtype=float)

        p_vec =  2 * beta_se_maf[indices_of_interest[j],2] * \
                    (1 - beta_se_maf[indices_of_interest[j],2]) * \
                    2*beta_se_maf[indices_of_interest,2] * \
                    (1 - beta_se_maf[indices_of_interest,2])
        w_vec = np.sum(reference_geno[:,j] ** 2) * w_sums
        covariances = reference_geno.transpose() @ reference_geno[:, j]

        new_b[j,:] = n_vec * np.sqrt(  p_vec / w_vec  ) * covariances

    new_b_inv = np.linalg.inv(new_b)
    beta = new_b_inv @ new_d @ beta_se_maf[indices_of_interest,0]


    return np.asarray([beta, np.sqrt(np.diag(new_b_inv))], dtype=float).transpose()