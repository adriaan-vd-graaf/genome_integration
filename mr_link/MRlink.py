import subprocess
import numpy as np
import scipy.stats
import argparse
import copy
import os
from genome_integration import resources
from genome_integration import simulate_mr
from genome_integration import causal_inference
from genome_integration.association import GeneticAssociation
from genome_integration import gene_regions
from genome_integration import utils

import time


def remove_plink_files(plink_file):

    for extension in ['.log', '.nosex', '.bed', '.bim', '.fam']:
        filename = f'{plink_file}{extension}'
        if os.path.exists(filename):
            os.remove(filename)


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
    sample_name_index_dict = {sample_names[x]: x for x in range(len(sample_names))}

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
                    phenotype_vector[sample_name_index_dict[sample_name]] = float(phenotype)
                    phenotype_available[sample_name_index_dict[sample_name]] = True
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
                        required=False,
                        default='NA',
                        help="ENSG id, so we can identify the region used and isolate "
                             "this from the reference and the outcome bed files."
                             "Please use the --region option to specify custom regions"
                        )

    parser.add_argument('--region',
                        type=str,
                        required=False,
                        default='NA',
                        help='A string that contains a genomic region. <chr>:<start_position>-<end_position>'
                             'Can be used as an alternative to --ensg_id.'
                             'An example would be to supply "2:1000-1100", which would indicate to take 1000 to 1100 on chromosome 2'
                             'this will later be padded with 1.5Mb of SNPs'
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

    mr_link_start_time = time.time()

    args = parser.parse_args()
    tmp_dir = args.temporary_location_prepend
    run_name = args.ensg_id if args.ensg_id != 'NA' else args.region
    big_outcome_bed_file = args.outcome_bed_file
    big_reference_bed_file = args.reference_bed_file

    tmp_subset_bed_outcome = f"{tmp_dir}{run_name}_subset_bed_outcome"
    tmp_subset_bed_reference = f"{tmp_dir}{run_name}_subset_bed_reference"
    ref_geno_for_cojo = f"{tmp_dir}{run_name}_subset_bed_reference_cojo"
    tmp_cojo = f"{tmp_dir}{run_name}_subset_bed_reference_tmp_cojo"
    outcome_pheno_file = args.outcome_phenotype_file

    gene_info = resources.read_gene_information()

    ## define the region that is used.
    if args.ensg_id == 'NA' and args.region == 'NA':
        raise ValueError('Please supply an argument to --ensg_id or to --region to determine the region that is analyzed.')
    elif args.ensg_id != 'NA' and args.region != 'NA':
        raise ValueError(
            'Please supply _ONE_ argument to _EITHER_ --ensg_id or to --region to determine the region that is analyzed.')
    elif args.ensg_id != 'NA' and args.region == 'NA':
        ##provided an ensembl ID.
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
    elif args.ensg_id == 'NA' and args.region != 'NA':
        # try to match region specification
        try:
            chr, _, start_end = args.region.partition(':')
            start, end = start_end.split('-')
            start = int(start)
            end = int(end)
        except Exception as x:
            raise ValueError(f'Could not parse {args.region} into the "<chr>:<start>-<end>" format.\n'
                             f'Exception message:{x}')

        region_of_interest = gene_regions.StartEndRegion([chr, start, end])
        region_size = region_of_interest.end - region_of_interest.start
        if region_of_interest.end - region_of_interest.start < 3000000:
            print(
                'WARNING: you specified a region that is smaller than 3Mb, for optimal results, please pad the edges of your region with a sufficient amount of SNPs.')
    else:
        raise ValueError('Programming logic error, please contact maintainer.')





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

    try:
        exposure_cojo = utils.do_gcta_cojo_on_genetic_associations(exposure_assocs,
                                                                   ref_geno_for_cojo,
                                                                   tmp_cojo,
                                                                   p_val_thresh=args.p_val_iv_selection_threshold,
                                                                   maf=0.01)

    except:
        raise ValueError(
            f"GCTA-COJO was unable to identify conditionally independent variants for {args.exposure_summary_statistics}")

    if len(exposure_cojo.ma_results.keys()) == 0:
        raise ValueError(
            f"GCTA-COJO found no independent variants for {args.exposure_summary_statistics} "
            f"at the {args.p_val_iv_selection_threshold:.2e} threshold"
        )

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
            exposure_cojo.ma_results[snp_name].beta *= -1
            exposure_cojo.ma_results[snp_name].z_score *= -1


    outcome_ld = simulate_mr.ld_from_geno_mat(outcome_geno)

    outcome_maf = np.apply_along_axis(simulate_mr.geno_frq, 1, outcome_geno)

    #runtime checks
    if reference_geno.shape[1] != outcome_geno.shape[1]:
        raise ValueError("Reference genotypes does not contain the same SNPs as outcome geno.")
    if outcome_plinkfile.bim_data.snp_names != reference_plinkfile.bim_data.snp_names:
        raise ValueError("The loci list for the outcome genotypes is not the same as the exposure genotypes")

    iv_selection_names = list(exposure_cojo.ma_results.keys())
    print(f'Debug: the names of the snps are {sorted(iv_selection_names)}')
    iv_selection = np.asarray([outcome_plinkfile.bim_data.snp_names.index(x) for x in iv_selection_names], dtype=int)

    scaled_reference_geno = np.apply_along_axis(simulate_mr.scale_geno_vec, 1, reference_geno)
    scaled_outcome_geno = np.apply_along_axis(simulate_mr.scale_geno_vec, 1, outcome_geno)

    iv_names = [reference_plinkfile.bim_data.snp_names[x] for x in iv_selection]
    beta_ses = np.asarray([[exposure_cojo.ma_results[x].beta, exposure_cojo.ma_results[x].se] for x in iv_names])

    mr_link_results = causal_inference.mr_link_ridge(
        scaled_outcome_geno,
        outcome_ld ** 2,
        beta_ses[:, 0],
        iv_selection,
        outcome_phenotypes,
        )

    print(f"Finished MR-link for {run_name} in {time.time() - mr_link_start_time} seconds.")

    with open(args.output_file, "w") as f:
        iv_summary_string = ','.join(
            [f'{exposure_cojo.ma_results[x].snp_name};'
             f'{exposure_cojo.ma_results[x].effect_allele};'
             f'{exposure_cojo.ma_results[x].beta:.5f};'
             f'{exposure_cojo.ma_results[x].se:.5f};'
             f'{exposure_cojo.ma_results[x].minor_allele_frequency:.3f}' for x in iv_names])
        f.write("ensembl_name\tmethod\tbeta\tse\tp_value\tn_ivs\tiv_summary\n")

        f.write(f"{run_name}\t"
                f"MR-link_uncalibrated\t"
                f"{mr_link_results[0]:.5f}\t"
                f"{mr_link_results[1]:.5f}\t"
                f"{mr_link_results[2]:.3e}\t"
                f"{len(exposure_cojo.ma_results)}\t"
                f"{iv_summary_string}\n")

    print(f"Uncalibrated MR-link results: beta: {mr_link_results[0]:.4f}, se: {mr_link_results[1]:.5f}, "
          f"p value: {mr_link_results[2]:.2e}")




