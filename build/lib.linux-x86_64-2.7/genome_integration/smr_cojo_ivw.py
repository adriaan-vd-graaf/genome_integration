import ivw_utils
import subprocess
import pickle
import scipy.stats
import isolate_snp_list_from_ma
import timeit
import time
import multiprocessing
import re
import file_utils
import numpy as np

from scipy.stats import pearsonr
from ma_utils import MaFile
from gcta_cojo_utils import CojoCmaFile
from variant_utils import BimFile

"""
These functions are being used to do inverse variance MR
of independent effects between gwas and eqtl results.
"""

__author__      = "Adriaan van der Graaf"
__copyright__   = "Copyright 2017, Adriaan van der Graaf"



class IVWFile:
    def __init__(self, file_name):
        lines = file_utils.read_newline_separated_file_into_list(file_name)[1:]
        self.ivw_lines = {}
        for line in lines:
            tmp = IVWLine(line)
            self.ivw_lines[tmp.probe_name] = tmp


    def find_probes_with_snp(self, snp_name):
        returnlist = set()
        for i in self.ivw_lines.keys():
            tmp = self.ivw_lines[i]
            if snp_name in tmp.snp_names:
                returnlist.add(tmp.probe_name)

        return returnlist


class IVWLine:
    def __init__(self, line):
        split = line.split("\t")
        self.probe_name = split[0]
        self.chromosome = split[1]
        self.n_snps = int(split[2])
        self.pearson_r = float(split[3])
        self.pearson_r_p = float(split[4])
        self.beta_ivw = float(split[5])
        self.se_ivw = float(split[6])
        self.p_ivw = float(split[7])
        self.snp_names = set(split[8].split(","))
        self.snp_smr_beta = split[9].split(",")
        self.snp_smr_se = split[10].split(",")


def do_gcta_cojo_slct(bfile_prepend, ma_file, out_prepend, p_val='1e-8', maf='0.01'):

    std_out = open(out_prepend + '.out', 'w')
    subprocess.run(['gcta64',
                    '--bfile', bfile_prepend,
                    '--cojo-file', ma_file,
                    '--cojo-slct',
                    '--out', out_prepend,
                    '--cojo-p', p_val,
                    '--maf', maf,
                    '--thread-num', '1'
                    ],
                   stdout=std_out,
                   check=True
                   )
    std_out.close()

    ##make sure the process is finished.
    regex = re.compile("^Analysis finished:.*")

    #make sure the log file is valid, assuming the other files are valid as well.
    for i in range(2):
        log_file = file_utils.read_newline_separated_file_into_list(out_prepend + '.out')
        # if len(log_file) == 0:
        #     time.sleep(1)
        #     continue
        if sum([regex.match(x) != None for x in log_file]) != 1:
            time.sleep(1)
        else:
            break
        if i == 9:
            raise IOError('Cojo analysis was finished, but did not find valid log file here:' + out_prepend + '.out' )

    return CojoCmaFile(out_prepend + ".jma.cojo", out_prepend)



def do_inverse_variance_weighting(probe_name,
                                  pval_thresh_all = '5e-2',
                                  pval_thresh_gwas = None,
                                  pval_thresh_eqtl = None,
                                  save_directory = "ivw_working_directory/",
                                  gwas_ma_loc_full ='Celiac_SMR_file.txt',
                                  ld_overlap_thres= 0.5,
                                  pickle_loc = None,
                                  full_plink_prepend = 'MAF0.01_04_INFO_ALL_including_chr10_only_overlap'
                                  ):
    """

    This Program will do inverse variance weighting for a list of genes.
    It will roughly work like this.

    1. Isolate snps from the bed file and 2 assoc statistic files (eqtl and gwas) for easy filtering
    2. Do the Cojo analysis on both association statistic files
    3. Determine the beta and se of SMR
    4. Do the inverse variance weighting.

    """

    if pval_thresh_eqtl == None:
        pval_thresh_eqtl = pval_thresh_all
    if pval_thresh_gwas == None:
        pval_thresh_gwas = pval_thresh_all


    print('Starting probe ' + probe_name, end='')

    # 1. To ensure analysis is relatively quick: isolate all the files which are in analysis
    # isolate the snps in the plink file, so that LD calculations are quicker
    # (LD calculations are done by GCTA COJO)

    probe_start = timeit.default_timer()

    probe_ma_loc = 'probes_for_besd/' + probe_name + '.ma'

    this_plink_file = save_directory + 'only_probe_' + probe_name + '_snps'
    snp_list_loc = save_directory + 'only_probe_' + probe_name + '_snps.txt'

    try:
        (ma_eqtl_data, bim_file) = isolate_snps_of_interest_make_bed(
            probe_ma_loc,
            probe_name,
            full_plink_prepend,
            snp_list_loc,
            this_plink_file
        )
    except subprocess.CalledProcessError:
        print("Could not find any snps overlapping for probe " + probe_name)
        return None

    ma_eqtl_data.add_bim_data(bim_file)

    # this returns nothing, but intersects both files, and creates the new file.
    gwas_ma_loc = save_directory + 'Celiac_SMR_file_only_' + probe_name

    ma_gwas_data = isolate_snp_list_from_ma.isolate_snps_from_list(
        snp_list_loc,
        gwas_ma_loc_full,
        gwas_ma_loc
    )

    ma_gwas_data.add_bim_data(bim_file)


    # now do the gcta cojo analysis.
    gcta_eqtl_cojo_out = save_directory + "gcta_initial_cojo_eqtl_" + probe_name
    try:
        eqtl_cojo_results = do_gcta_cojo_slct(this_plink_file,
                                              probe_ma_loc,
                                              gcta_eqtl_cojo_out,
                                              p_val=pval_thresh_eqtl
                                              )

        #not needed anymore, but still implemented.
        #cojo_ld_data = CojoLdrFile(gcta_eqtl_cojo_out + ".ldr.cojo", gcta_eqtl_cojo_out)


    except FileNotFoundError:
        print("\rCojo did not yield any results for probe " + probe_name)
        print("\tcontinueing with the next probe")
        return None

    # and do it for the GWAS data
    gcta_gwas_cojo_out = save_directory + "gcta_initial_cojo_gwas_" + probe_name

    try:
        gwas_cojo_results = do_gcta_cojo_slct(this_plink_file,
                                              gwas_ma_loc,
                                              gcta_gwas_cojo_out,
                                              p_val=pval_thresh_gwas
                                              )
    except FileNotFoundError:
        print("Warning: Cojo did not yield any results for the file " + gwas_ma_loc)
        print("\tcontinueing with the next probe")
        return None

    snps = bim_file.isolate_LD_similar_snps(gwas_cojo_results.snps_with_data(),
                                            eqtl_cojo_results.snps_with_data(),
                                            ld_overlap_thres
                                            )

    ivw_result = ivw_utils.IVWResult()
    # initialize the list for all snps.

    outcome_z = [0.0] * len(snps)
    exposure_z = [0.0] * len(snps)

    for i in range(len(snps.keys())):
        snp = list(snps.keys())[i]

        gwas_snp = snp
        eqtl_snp, best_ld = snps[snp]

        eqtl_cojo_results.ma_results[eqtl_snp].correct_score_based_on_ld(best_ld)

        exposure_z[i] = gwas_cojo_results.ma_results[gwas_snp].get_z_score()
        outcome_z[i] = eqtl_cojo_results.ma_results[eqtl_snp].get_z_score()

        smr_result = estimate_beta_se_smr(
                                            gwas_cojo_results.ma_results[gwas_snp],
                                            eqtl_cojo_results.ma_results[eqtl_snp]
                                        )

        # bit of a hack, but if the sign is different in the LD,
        # then the sign of the beta_smr should also be different.

        smr_result = (smr_result[0] * np.sign(best_ld), smr_result[1], smr_result[2])

        try:
            pos_tmp = bim_file.bim_results[snp].position()
            chr_tmp = bim_file.bim_results[snp].chromosome()
        except:
            pos_tmp = np.NaN
            chr_tmp = np.NaN

        ivw_result.add_estimate(smr_result,
                                gwas_snp,
                                pos_tmp,
                                chr_tmp,
                                eqtl_snp
                            )

    #determine correlation of the z_score estimates
    if len(outcome_z) > 1:
        z_score_correlation, z_score_correlation_p = pearsonr(outcome_z, exposure_z)
    else:
        #make sure a warning is not posted.
        z_score_correlation, z_score_correlation_p = (np.nan, np.nan)


    # now do the inverse variance weighting.
    try:
        ivw_result.do_ivw_estimation()
    except:
        print("\rWarning: Zero variance of estimators found for probe " + probe_name)
        print("\tContinueing with the next probe")
        return None

    try:
        (beta_ivw, se_ivw, p_value_ivw) = ivw_result.get_ivw_estimates()
    except RuntimeError:
        print("\rNo beta smr estimates provided in IVW estimation" + probe_name)
        return None


    # write away results
    probe_end = timeit.default_timer()
    print("\rFinished probe " + probe_name
          + ' in {:4.2f}'.format(probe_end - probe_start)
          + ' seconds')

    #save results to a pickle, if we need more plotting
    for_pickle = {"ivw_data": ivw_result,
                  "ma_eqtl": ma_eqtl_data,
                  "ma_gwas": ma_gwas_data,
                  "cojo_eqtl": eqtl_cojo_results,
                  "gwas_cojo": gwas_cojo_results,
                  "bim_data": bim_file
                  }

    if pickle_loc == None:
        pickle_loc = save_directory + probe_name + '_pickle.p'

    pickle.dump(for_pickle, open(pickle_loc, 'wb'))


    #now remove all files.
    subprocess.run(['rm',
                    save_directory + 'gcta_initial_cojo_eqtl_' + probe_name + '.jma.cojo',
                    save_directory + 'gcta_initial_cojo_eqtl_' + probe_name + '.ldr.cojo',
                    save_directory + 'gcta_initial_cojo_eqtl_' + probe_name + '.out',
                    save_directory + 'gcta_initial_cojo_gwas_' + probe_name + '.jma.cojo',
                    save_directory + 'gcta_initial_cojo_gwas_' + probe_name + '.ldr.cojo',
                    save_directory + 'gcta_initial_cojo_gwas_' + probe_name + '.out',
                    save_directory + 'only_probe_'+ probe_name + '_snps.txt',
                    save_directory + 'only_probe_'+ probe_name + '_snps.bed',
                    save_directory + 'only_probe_'+ probe_name + '_snps.bim',
                    save_directory + 'only_probe_'+ probe_name + '_snps.fam',
                    save_directory + 'only_probe_'+ probe_name + '_snps.log',
                    save_directory + 'Celiac_SMR_file_only_' + probe_name
                    ])

    return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                    probe_name,
                                    gwas_cojo_results.ma_results[list(snps.keys())[0]].chr,
                                    len(snps),
                                    z_score_correlation,
                                    z_score_correlation_p,
                                    beta_ivw,
                                    se_ivw,
                                    p_value_ivw,
                                    ','.join(snps),
                                    ','.join([str(x[0]) for x in ivw_result.estimation_data]),
                                    ','.join([str(x[1]) for x in ivw_result.estimation_data])
                                    )


def main():
    probe_names = file_utils.read_newline_separated_file_into_list("interesting_probes.txt")

    # make sure there is a directory to which we write
    save_directory = "ivw_working_directory/"
    subprocess.run(['mkdir', '-p', save_directory], check=True)
    # the GWAS data

    thread_pool = multiprocessing.Pool(4)
    output_list = thread_pool.map(do_inverse_variance_weighting, probe_names)

    write_file = open("ivw_results.txt", 'w')
    write_file.write(
        "probe_name\tchromosome\tn_snps\tpearson_r\tpearson_p\tbeta_ivw\tse_ivw\tp_value_ivw\tsnp_names\tsnp_smr_betas\tsnp_smr_se\n")

    for i in output_list:
        if i != None:
            write_file.write(i)

    write_file.close()

if __name__ == "__main__":
    main()