import subprocess
import numpy as np
import sqlite3
import scipy.stats
import argparse
import copy
from genome_integration import resources
from genome_integration import simulate_mr
from genome_integration import causal_inference
from genome_integration.association import GeneticAssociation
from genome_integration import gene_regions
from genome_integration import utils

import time


def remove_plink_files(plink_file):
    subprocess.run(
        ["rm {}.*".format(plink_file)],
        shell=True,
        check=True,
        # stdout=subprocess.DEVNULL,
        # stderr=subprocess.DEVNULL:
    )


def check_db(gene_name):
    db = sqlite3.connect("gtex_analysis/dbs/" + gene_name + "_ldl_cholesterol_results.db", timeout=600)

    try:
        db.execute('''CREATE TABLE ldl_cholesterol_results_003
                        (tissue text,
                         gene_name text, 
                         method text, 
                         selection_method text, 
                         p_val_thresh real, 
                         estimated_beta real, 
                         estimated_se real, 
                         p_val_of_estimate real,
                         snps_identified real
                         )''')

    except Exception as x:
        print(x)

    db.commit()
    db.close()


def save_result_to_db(beta_se_n_tuple, method, selection_method, gene_name, threshold, tissue):
    database = sqlite3.connect("gtex_analysis/dbs/" + gene_name + "_ldl_cholesterol_results.db", timeout=600)

    database.execute("""INSERT  INTO ldl_cholesterol_results_003 VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                     (tissue,
                      gene_name,
                      method,
                      selection_method,
                      threshold,
                      beta_se_n_tuple[0],
                      beta_se_n_tuple[1],
                      beta_se_n_tuple[3],
                      beta_se_n_tuple[2],
                      )
                     )

    database.commit()
    database.close()

def read_outcome_genotypes_and_phenotypes(big_bed, tmp_loc, region,phenotype_file=None, variants_to_keep=None):
    bed_file = None
    genotypes = None
    phenotype_vector = None

    # here we start making some files, but we WILL clean up after ourselves using the remove plink file function.

    utils.read_region_from_plink(big_bed, tmp_loc, region, variants_to_keep)

    bed_file = utils.PlinkFile(tmp_loc)
    genotypes = bed_file.read_bed_file_into_numpy_array()

    sample_names = [x for x in bed_file.fam_data.fam_samples]

    # geno matrix is read, now read phenotypes, and
    # then prune the genotypes for the individuals without phenotype information

    phenotype_available = np.zeros(len(sample_names), dtype=bool)
    phenotype_vector = np.zeros(len(sample_names), dtype=float)

    if phenotype_file is not None:
        n = 0
        correct_header="FID\tIID\tPHENO\n"
        with open(phenotype_file) as f:
            header= f.readline()
            if header != correct_header:
                raise ValueError("Outcome phenotype file should have the correct tab separated newline ended header:",
                                 header)
            for line in f:
                fid, iid, phenotype = line.split()
                sample_name = f"{fid}~__~{iid}"
                try:
                    phenotype_vector[sample_names.index(sample_name)] = float(phenotype)
                    phenotype_available[sample_names.index(sample_name)] = True
                except:
                    n += 1
        if n > 0:
            print('could not find', n, "iids in genotypes")

        # now prune the genotypes, and phenotypes.
        genotypes = genotypes[phenotype_available,:]
        sample_names_to_remove = np.asarray(sample_names)[~phenotype_available]
        for sample_to_remove in sample_names_to_remove:
            del bed_file.fam_data.fam_samples[sample_to_remove]
            bed_file.fam_data.sample_names.remove(sample_to_remove)

        phenotype_vector = phenotype_vector[phenotype_available]

    remove_plink_files(tmp_loc)

    return bed_file, genotypes, phenotype_vector


def read_summary_statistic_file(file_name):
    assocs = {}
    correct_header = "CHR\tPOS\tNAME\tREF_ALLELE\tEFFECT_ALLELE\tBETA\tSE\tMAF\tN_OBS\n"
    with open(file_name, "r") as f:
        header= f.readline()
        if header != correct_header:
            raise ValueError(f"Summary statistic file should be the same as the following\n"
                             f"{correct_header}"
                             f"all fields tab separated and newline ended.")
        for line in f:
            split = line.split()
            tmp_assoc = GeneticAssociation(dependent_name=args.ensg_id,
                                             explanatory_name=split[2],
                                             n_observations=split[8],
                                             beta=split[5],
                                             se=split[6],
                                             r_squared=None,
                                             chromosome=split[0],
                                             position=split[1],
                                             major_allele=split[3],
                                             minor_allele=split[4],
                                             minor_allele_frequency=split[7],
                                             reference_allele=None,
                                             effect_allele=None
                                             )
            tmp_assoc.set_p_val(scipy.stats.t.sf(np.abs(tmp_assoc.z_score), tmp_assoc.n_observations-2) * 2)

            assocs[split[2]] = tmp_assoc


    return assocs




if __name__ == '__main__':

    total_start_time = time.time()
    parser = argparse.ArgumentParser()

    parser.add_argument("--outcome_bed_file",
                        type=str,
                        required=True,
                        help="plink bed file of the outcome, in b37, this will be pruned for the region of the ensg_id")

    parser.add_argument("--outcome_phenotype_file",
                        type=str,
                        required=True,
                        help="phenotype file of the outcome")

    parser.add_argument("--reference_bed_file",
                        type=str,
                        required=True,
                        help="plink bed file of the reference cohort, in b37, "
                             "can be the same as the outcome bed file, but results are untested")

    parser.add_argument("--exposure_summary_statistics",
                        type=str,
                        required=True,
                        help="eQTL summary statistics file of the exposure summary statistics"
                        )
    parser.add_argument("--ensg_id",
                        type=str,
                        required=True,
                        help="ENSG id, so we can identify the region used and isolate "
                             "this from the reference and the outcome bed files."
                        )

    parser.add_argument("--temporary_location_prepend",
                        type=str,
                        default="temp_mr_link_files_",
                        help="location where temporary files can be stored.")

    parser.add_argument("--p_val_iv_selection_threshold",
                        default=5e-8,
                        type=float,
                        help="p value selection theshold for GCTA COJO"
                        )

    parser.add_argument("--output_file",
                        type=str,
                        help="the file name where the MR-link results will be output."
                        )

    parser.add_argument("--permute", action='store_true',
                        help="Use this flag if you want permutations to occur (will add substantial runtime)")

    mr_link_start_time = time.time()

    args = parser.parse_args()
    tmp_dir = args.temporary_location_prepend
    ensg_name = args.ensg_id
    big_outcome_bed_file = args.outcome_bed_file
    big_reference_bed_file = args.reference_bed_file

    tmp_subset_bed_outcome = f"{tmp_dir}{ensg_name}_subset_bed_outcome"
    tmp_subset_bed_reference = f"{tmp_dir}{ensg_name}_subset_bed_reference"
    ref_geno_for_cojo = f"{tmp_dir}{ensg_name}_subset_bed_reference_cojo"
    tmp_cojo = f"{tmp_dir}{ensg_name}_subset_bed_reference_tmp_cojo"
    outcome_pheno_file = args.outcome_phenotype_file

    gene_info = resources.read_gene_information()

    if args.ensg_id == "simulated_run" or args.ensg_id == 'ENSG00000000000':
        ensg_info = gene_regions.StartEndRegion(["2", 100000000, 105000000])

        region_of_interest = copy.deepcopy(ensg_info)
    else:
        # start end region.
        ensg_info = gene_info.str_to_full(args.ensg_id)

        try:
            region_of_interest = gene_regions.StartEndRegion([ensg_info.chromosome,
                                                              ensg_info.start - 1500000,
                                                              ensg_info.end + 1500000])
        except:
            region_of_interest = gene_regions.StartEndRegion([ensg_info.chromosome,
                                                              0,
                                                              ensg_info.end + 1500000])

    #read in exposure summary statistics
    exposure_assocs = read_summary_statistic_file(args.exposure_summary_statistics)


    #read in the outcome genotypes.
    outcome_plinkfile, outcome_geno, outcome_phenotypes = read_outcome_genotypes_and_phenotypes(big_outcome_bed_file,
                                                                                               tmp_subset_bed_outcome,
                                                                                               region_of_interest,
                                                                                               outcome_pheno_file,
                                                                                               exposure_assocs.keys()
                                                                                               )

    #read in the reference genotypes.
    reference_plinkfile, _, _ = read_outcome_genotypes_and_phenotypes(big_reference_bed_file,
                                                                        tmp_subset_bed_reference,
                                                                        region_of_interest,
                                                                        phenotype_file=None,
                                                                        variants_to_keep=exposure_assocs.keys()
                                                                        )

    #harmonization of the genotypes.
    reference_plinkfile, reference_geno = outcome_plinkfile.harmonize_genotypes(reference_plinkfile)

    #run GCTA COJO. there is duplication here, but handling files and their deletion can be hard.
    utils.read_region_from_plink(big_reference_bed_file,
                                 ref_geno_for_cojo,
                                 region_of_interest,
                                 variants=exposure_assocs.keys())

    exposure_cojo = utils.do_gcta_cojo_on_genetic_associations(exposure_assocs,
                                                               ref_geno_for_cojo,
                                                               tmp_cojo,
                                                               p_val_thresh=args.p_val_iv_selection_threshold,
                                                               maf=0.01)

    #remove the isolated bed, as it's not being removed yet.
    subprocess.run([f"rm {ref_geno_for_cojo}.*"], shell=True, check=True)

    #Harmonize COJO alleles to the outcome.
    for snp_name in exposure_cojo.ma_results.keys():
        if (exposure_cojo.ma_results[snp_name].minor_allele !=
            outcome_plinkfile.bim_data.bim_results[snp_name].minor_allele) and (
                exposure_cojo.ma_results[snp_name].minor_allele !=
                outcome_plinkfile.bim_data.bim_results[snp_name].major_allele):
            raise ValueError("Alleles did not match between ")

        elif (exposure_cojo.ma_results[snp_name].minor_allele !=
            reference_plinkfile.bim_data.bim_results[snp_name].minor_allele):
            #flip the alleles and effect size:
            old_major = exposure_cojo.ma_results[snp_name].major_allele
            exposure_cojo.ma_results[snp_name].major_allele = exposure_cojo.ma_results[snp_name].minor_allele
            exposure_cojo.ma_results[snp_name].minor_allele = old_major
            exposure_cojo.beta *= -1
            exposure_cojo.z_score *= -1


    outcome_ld = simulate_mr.ld_from_geno_mat(outcome_geno)

    outcome_maf = np.apply_along_axis(simulate_mr.geno_frq, 1, outcome_geno)

    #runtime checks
    if reference_geno.shape[1] != outcome_geno.shape[1]:
        raise ValueError("Reference genotypes does not contain the same SNPs as outcome geno.")
    if outcome_plinkfile.bim_data.snp_names != reference_plinkfile.bim_data.snp_names:
        raise ValueError("The loci list for the outcome genotypes is not the same as the exposure genotypes")

    iv_selection_names = list(exposure_cojo.ma_results.keys())
    iv_selection = np.asarray([outcome_plinkfile.bim_data.snp_names.index(x) for x in iv_selection_names], dtype=int)


    scaled_reference_geno = np.apply_along_axis(simulate_mr.scale_geno_vec, 1, reference_geno)
    scaled_outcome_geno = np.apply_along_axis(simulate_mr.scale_geno_vec, 1, outcome_geno)

    iv_names = [reference_plinkfile.bim_data.snp_names[x] for x in iv_selection]
    beta_ses = np.asarray([[exposure_cojo.ma_results[x].beta, exposure_cojo.ma_results[x].se] for x in iv_names])
    #
    # beta_ses = simulate_mr.do_gcta_cojo_conditional(scaled_reference_geno,
    #                                              exposure_assocs,
    #                                              iv_selection,
    #                                              outcome_plinkfile.bim_data.snp_names
    #                                              )

    mr_link_results = causal_inference.mr_link_ridge(
        scaled_outcome_geno,
        outcome_ld ** 2,
        beta_ses[:, 0],
        iv_selection,
        outcome_phenotypes,
        )

    permuted_p = np.nan
    print(f"Finished MR-link for {ensg_name} in {time.time() - mr_link_start_time} seconds.")
    #permutation scheme.
    if args.permute:
        print("Starting permutations")
        permuted_p_values = []
        n_permutations = 1000
        for i in range(n_permutations):
            permuted_p_values.append(causal_inference.mr_link_ridge(
                scaled_reference_geno,
                outcome_ld ** 2,
                beta_ses[:, 0],
                iv_selection,
                np.random.permutation(outcome_phenotypes),
            )[2])

        permuted_p = np.sum(permuted_p_values < mr_link_results[2]) / n_permutations

        print(f"Finished MR-link and permutations for {ensg_name} in {time.time() - mr_link_start_time} seconds.")


    with open(args.output_file, "w") as f:
        iv_summary_string = ','.join(
            [f'{exposure_cojo.ma_results[x].snp_name};'
             f'{exposure_cojo.ma_results[x].effect_allele};'
             f'{exposure_cojo.ma_results[x].beta:.5f};'
             f'{exposure_cojo.ma_results[x].se:.5f};'
             f'{exposure_cojo.ma_results[x].minor_allele_frequency:.3f}' for x in iv_names])
        f.write("ensembl_name\tmethod\tbeta\tse\tp_value\tn_ivs\tiv_summary\n")

        f.write(f"{ensg_name}\t"
                f"MR-link_uncalibrated\t"
                f"{mr_link_results[0]:.5f}\t"
                f"{mr_link_results[1]:.5f}\t"
                f"{mr_link_results[2]:.3e}\t"
                f"{len(exposure_cojo.ma_results)}\t"
                f"{iv_summary_string}\n")

        if args.permute:
            f.write(f"{ensg_name}\t"
                    f"MR-link_permuted_iids\t"
                    f"{mr_link_results[0]:.5f}\t"
                    f"{mr_link_results[1]:.5f}\t"
                    f"{permuted_p:.3e}\t"
                    f"{len(exposure_cojo.ma_results)}\t"
                    f"{iv_summary_string}\n")


    print(f"Uncalibrated MR-link results: beta: {mr_link_results[0]:.4f}, se: {mr_link_results[1]:.5f}, p value: {mr_link_results[2]:.2e}")




