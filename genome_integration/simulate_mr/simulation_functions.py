import copy
from scipy.stats import binom
import numpy as np
import statsmodels.api as sm

from . import geno_functions


def simulate_phenotypes_binary_outcome(
        exposure_1_causal, exposure_2_causal,
        exposure_1_n_causal, exposure_2_n_causal,
        exposure_geno, exposure_ld,
        outcome_geno,
        exposure_ld_r_sq = None,
        overlapping_causal_snps = 0,
        error_sd = 1.0,
        confounder_sd = 0.5,
        inside_phi = 0.0,
        directional_pleiotropy = True,
        known_exposure_lower_ld_bound=0.0,
        proportion_cases=0.5
        ):
    """

    This function simulates two binary outcome phenotypes in cohorts which are under the same genetic control

    Step by step:

    1.
        exposure 1 SNPs are chosen, and they are chosen not to be in higher ld than 0.95 R ^ 2

    2.
        exposure 2 SNPs are chosen in two ways. `overlapping_causal_snps` are chosen without replacement from the
        exposure 1 SNPs, the rest are chosen from SNPs that are in LD with exposure 1.

    3.
        The confounder is simulated, based on some error and if `inside_phi` is non-zero, it is given a genetic effect.

    4.
        The effect sizes of exposure 1 and of exposure 2 eQTLs are chosen.

    5.
        The phenotypes are simulated based on the function simulate_phenotypes_extended,
        in the exposure cohort, and the outcome cohort independently.

    6.
        The outcome phenotypes are transformed using a logistic function. 
        The they are converted into case control variables following the binomial function with the logistically transformed variables as probabilities.
    
    7.
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
    The ld R matrix for the exposure genotypes

    :param outcome_geno:
    The genotype matrix for the outcome encoded

    :param outcome_ld:
    The ld R matrix for the outcome genotypes.

    :param overlapping_causal_snps: int
    An integer number representing the number of causal SNPs for exposure 2 that overlap exposure 1.

    :param error_sd: float
    The scale (standard deviation) parameter of the error that is added to a phenotype

    :param inside phi: float
    The phi parameter (from the egger paper notation), that makes the confounder dependent of the genotype.
    if non zero, a genetic effect with effect size (0,2) is added to the confounder, using the exposure 2 causal SNPs

    :param directional_pleiotropy: bool
    if there should be directional pleiotropy or not.

    :param proportion_cases: float
        The proportion of cases and for the outcome phenotypes.


    :return:
    Returns a tuple of the following numpy arrays:

    1 by n vector of simulated phenotypes for the exposure,
    1 by n vector simulated binary phenotypes for the outcome, 1 = case, 0 = control
    1 by n causal vector of the indices of the causal snps of exposure 1.
    1 by n causal vector of the effect sizes of the causal snps of exposure 1.
    1 by n causal vector of the indices of the causal snps of exposure 2.
    1 by n causal vector of the effect sizes of the causal snps of exposure 2.

    """

    exposure_phenotype, outcome_phenotype, \
        exposure_1_causal_snps, exposure_1_betas, \
        exposure_2_causal_snps, exposure_2_betas =  simulate_phenotypes(
                                                        exposure_1_causal, exposure_2_causal,
                                                        exposure_1_n_causal, exposure_2_n_causal,
                                                        exposure_geno, exposure_ld,
                                                        outcome_geno,
                                                        exposure_ld_r_sq,
                                                        overlapping_causal_snps,
                                                        error_sd,
                                                        confounder_sd,
                                                        inside_phi,
                                                        directional_pleiotropy,
                                                        known_exposure_lower_ld_bound
                                                        )

    outcome_phenotype_logistic = np.asarray(
            [1 / (1 + np.exp(-x)) for x in outcome_phenotype], 
            dtype=float)
    
    outcome_phenotype = binom.rvs(1, outcome_phenotype_logistic)

    return exposure_phenotype, outcome_phenotype, \
           exposure_1_causal_snps, exposure_1_betas, \
           exposure_2_causal_snps, exposure_2_betas




def simulate_phenotypes(exposure_1_causal, exposure_2_causal,
                        exposure_1_n_causal, exposure_2_n_causal,
                        exposure_geno, exposure_ld,
                        outcome_geno,
                        exposure_ld_r_sq = None,
                        overlapping_causal_snps = 0,
                        error_sd = 1.0,
                        confounder_sd = 0.5,
                        inside_phi = 0.0,
                        directional_pleiotropy = True,
                        known_exposure_lower_ld_bound=0.0
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
        The confounder is simulated, based on some error and if `inside_phi` is non-zero, it is given a genetic effect.

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
    The ld R matrix for the exposure genotypes -- NB. NOT squared correlation

    :param outcome_geno:
    The genotype matrix for the outcome

    :param outcome_ld:
    The ld R matrix for the outcome genotypes. currently not used.

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


    if exposure_ld is not None:
        exposure_ld = exposure_ld ** 2
    elif exposure_ld_r_sq is not None and exposure_ld is None:
        exposure_ld = exposure_ld_r_sq
    else:
        raise ValueError("Programmer Error")

    """
    Step 1.
    """

    #make sure the exposure Causal SNPs are not in > upper_ld_bound R^2 with another.
    exposure_1_causal_snps = np.zeros([exposure_1_n_causal], dtype=int)
    i=0
    exposure_1_causal_snps[i] = np.random.choice(exposure_geno.shape[1], 1, replace=False)
    i += 1

    #i indicates how many snps have been chosen
    while i < exposure_1_n_causal:
        # choosing SNPS that are not in full linkage (R**2 > 0.95) with one another, but do share some linkage.
        lower_than_upper = np.sum(exposure_ld[exposure_1_causal_snps[:i], :] < upper_ld_bound, axis=0) == i

        if known_exposure_lower_ld_bound >= 0.0: #this part tries to find causal exposure variants in high LD.
            higher_than_lower = np.sum(exposure_ld[exposure_1_causal_snps[:i], :] >= known_exposure_lower_ld_bound, axis=0) >= 1
            if sum(higher_than_lower) == 0:
                raise ValueError("Could not find known exposure variants within the lower LD bound")
            if sum(higher_than_lower & lower_than_upper) == 0:
                raise ValueError("Could not find known exposure variants within the lower AND higher LD bound")
            choice = np.random.choice(np.where(higher_than_lower & lower_than_upper)[0], 1, replace=False)
        else:
            choice = np.random.choice(np.where(lower_than_upper)[0], 1, replace=False)

        exposure_1_causal_snps[i] = choice
        i += 1

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

        exposure_1_causal_snps_non_overlapping = np.asarray(
            [x for x in exposure_1_causal_snps if x not in exposure_2_overlapping_snps],
            dtype=int)

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

    confounder_exposure_cohort = np.random.normal(0, confounder_sd, exposure_geno.shape[0])
    confounder_outcome_cohort = np.random.normal(0, confounder_sd, outcome_geno.shape[0])

    if inside_phi != 0.0:
        phi_values = np.random.uniform(0, inside_phi, (1, len(exposure_2_causal_snps)))
        confounder_exposure_cohort += (exposure_geno[:,exposure_2_causal_snps] @ phi_values).reshape(exposure_geno.shape[0])
        confounder_outcome_cohort += (outcome_geno[:,exposure_2_causal_snps] @phi_values).reshape(outcome_geno.shape[0])


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
    exposure 1: y_{e1} = X_{1} cdot b_{e1} + C + e
    exposure 2: y_{e2} = X_{2} cdot b_{e2} + C + e

    outcome:  y_{o} =  y_{e1} b_{1} + X_{1} cdot b_{p,1} + X_{2} cdot b_{p,1} + y_{e2} b_{2} + C + e

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

    scaled_geno = np.apply_along_axis(geno_functions.scale_geno_vec, 0, geno)

    # exposure 1: y_{e1} = X_{1}b_{e1} + C + e
    exposure_1_phenotype = scaled_geno[:,exposure_1_causal_snps] @ exposure_1_betas

    exposure_1_phenotype += confounder
    exposure_1_phenotype += np.random.normal(0, error_scale, scaled_geno.shape[0])  # add some measurement error

    # exposure 2: y_{e2} = X_{2}b_{e2} + C + e
    exposure_2_phenotype = scaled_geno[:,exposure_2_causal_snps] @ exposure_2_betas
    exposure_2_phenotype += confounder
    exposure_2_phenotype += np.random.normal(0, error_scale, scaled_geno.shape[0])  # add some measurement error

    # y_{o} = y_{e1} b_{1} + X_{1}b_{p,1} + X_{2}_b_{p,1} + y_{e2} b_{2} + C + e
    outcome_phenotype = exposure_1_phenotype * exposure_1_causal_effect
    outcome_phenotype += exposure_2_phenotype * exposure_2_causal_effect

    outcome_phenotype += confounder
    outcome_phenotype += np.random.normal(0, error_scale, scaled_geno.shape[0])  # add some measurement error

    """
    Phenotypes are now simulated  
    """

    return exposure_1_phenotype, exposure_2_phenotype, outcome_phenotype




