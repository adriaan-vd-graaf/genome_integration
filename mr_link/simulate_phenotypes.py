import pickle
import os
import argparse
import numpy as np
import scipy.linalg
import scipy.stats

from genome_integration import simulate_mr
from genome_integration import utils
from genome_integration.association import GeneticAssociation

def turn_assocs_into_genetic_associations(assocs, ordered_loci, allele_frequency, sample_size):
    z_scores = assocs[:,0] / assocs[:,1]
    p_values = scipy.stats.norm.sf(np.abs(z_scores)) * 2

    assocs = {ordered_loci[i].name:
            GeneticAssociation(dependent_name="simulation",
                explanatory_name=ordered_loci[i].name,
                n_observations = sample_size,
                beta=assocs[i,0],
                se=assocs[i,1],
                r_squared = None,
                chromosome = ordered_loci[i].chromosome,
                position = ordered_loci[i].bp_position,
                major_allele = ordered_loci[i].allele2,
                minor_allele = ordered_loci[i].allele1,
                minor_allele_frequency = allele_frequency[i],
                reference_allele = None,
                effect_allele = None
                )
          for i in range(len(assocs))
          }
    [assocs[ordered_loci[i].name].set_wald_p_value(p_values[i]) for i in range(len(assocs))]
    return assocs


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("--bed_cohort_1",
                        required=True,
                        help="The bed file for cohort 1 (from which the exposure is calculated and saved)",

                        )
    parser.add_argument("--bed_cohort_2",
                        required=True,
                        help="The bed file for cohort 2 (from which the outcome is calculated and saved)",
                        )


    parser.add_argument("--out_prepend",
                        default="temp_simulated_file_",
                        help="The prepend for where to place the file")

    parser.add_argument("--save_as",
                        default="numpy",
                        help="How to save the results. "
                             "'numpy' will be more efficient and saves all data to a single file."
                             "'text' will be human readable but save to eight separate files per simulation." 
                             "The latter can blow up the number of files in the file system.",
                        choices=['numpy', 'text'],
                        # required=True
                        )

    parser.add_argument("--exposure_1_causal",
                        type=float,
                        default=0.0,
                        help="causal effect of exposure 1")

    parser.add_argument("--exposure_2_causal",
                        type=float,
                        default=0.4,
                        help="causal effect of the pleiotropic exposure 2")

    parser.add_argument("--n_causal_exposure_1",
                        type=int,
                        default=5,
                        help="number of causal SNPs for exposure 1")

    parser.add_argument("--n_causal_exposure_2",
                        type=int,
                        default=5,
                        help="number of causal SNPs for the pleiotropic exposure 2")

    parser.add_argument("--overlapping_causal",
                        type=int,
                        default=0,
                        help="number of causal SNPs that overlap between exposure 1 and exposure 2")

    parser.add_argument("--directional_pleiotropy", type=int, default=1,
                        help="If the effects of the second exposure should be directional")

    parser.add_argument("--phenotypes_to_simulate", type=int, default=5,
                        help= "How many phenotypes are saved could take a long time. But results are saved intermittendly")

    args = parser.parse_args()


    if args.save_as == "text":
        raise NotImplementedError("Text file saving is still a todo")


    sim_scenario = "ex1-n-{}_ex2-n-{}_ex1-b-{}_ex2-b-{}_overl-{}_dir_pleio-{}".format(
        args.n_causal_exposure_1,
        args.n_causal_exposure_2,
        args.exposure_1_causal,
        args.exposure_2_causal,
        args.overlapping_causal,
        args.directional_pleiotropy,
    )


    exposure_plink_file = utils.PlinkFile(args.bed_cohort_1)
    geno_exposure = exposure_plink_file.read_bed_file_into_numpy_array()

    outcome_plink_file = utils.PlinkFile(args.bed_cohort_1)
    geno_outcome = outcome_plink_file.read_bed_file_into_numpy_array()

    exposure_ld = simulate_mr.ld_from_geno_mat(geno_exposure)
    exposure_maf = np.apply_along_axis(simulate_mr.geno_frq, 0, geno_exposure)
    outcome_maf = np.apply_along_axis(simulate_mr.geno_frq, 0, geno_outcome)

    print("Reading genotype files finished")

    i = 0
    while i < args.phenotypes_to_simulate:

        file_name = f"{args.out_prepend}{sim_scenario}_run_{i}"

        try:

            exposure_phenotype, outcome_phenotype, \
            exposure_1_causal_snps, exposure_1_betas, \
            exposure_2_causal_snps, exposure_2_betas = \
                simulate_mr.simulate_phenotypes(args.exposure_1_causal, args.exposure_2_causal,
                                                         args.n_causal_exposure_1, args.n_causal_exposure_2,
                                                         geno_exposure, exposure_ld,
                                                         geno_outcome,
                                                         overlapping_causal_snps=args.overlapping_causal,
                                                         inside_phi=0.0,
                                                         directional_pleiotropy=args.directional_pleiotropy)


            #runtime check that the simulate phenotype function performs as expected
            iteration = 0
            while len(set(exposure_1_causal_snps) & set(exposure_2_causal_snps)) != args.overlapping_causal:
                exposure_phenotype, outcome_phenotype, \
                exposure_1_causal_snps, exposure_1_betas, \
                exposure_2_causal_snps, exposure_2_betas = \
                    simulate_mr.simulate_phenotypes(args.exposure_1_causal, args.exposure_2_causal,
                                                             args.n_causal_exposure_1, args.n_causal_exposure_2,
                                                             geno_exposure, exposure_ld,
                                                             geno_outcome,
                                                             overlapping_causal_snps=args.overlapping_causal,
                                                             directional_pleiotropy=args.directional_pleiotropy)
                iteration +=1
                if iteration > 100:
                    raise RuntimeError("tried 100 times to simulate non overlapping causal snps, failed all times.")

        except Exception as x:
            print("Failed simulating genotypes", i, "with error",  x )
            continue

        max_ld_of_causal = np.max(np.tril(exposure_ld[:,exposure_1_causal_snps][exposure_1_causal_snps,:] **2 ,-1))

        print("Simulated", i, "max ld between exposure 1 causal variants is: {:.3f}".format(max_ld_of_causal))


        exposure_sum_stats = np.apply_along_axis(simulate_mr.do_gwas_on_scaled_variants,
                                                 0, geno_exposure, dependent=exposure_phenotype).T

        outcome_sum_stats = np.apply_along_axis(simulate_mr.do_gwas_on_scaled_variants,
                                            0, geno_outcome, dependent=outcome_phenotype).T


        np.savez(file_name,
                 exposure_phenotype=exposure_phenotype,
                 outcome_phenotype=outcome_phenotype,
                 exposure_1_causal_snps=exposure_1_causal_snps,
                 exposure_1_betas=exposure_1_betas,
                 exposure_2_causal_snps=exposure_2_causal_snps,
                 exposure_2_betas=exposure_2_betas,
                 exposure_sum_stats=exposure_sum_stats,
                 outcome_sum_stats=outcome_sum_stats,
                )

        i += 1