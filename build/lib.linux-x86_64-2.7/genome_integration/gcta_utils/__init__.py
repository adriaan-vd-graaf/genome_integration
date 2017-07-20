"""
These classes are used to store and deal with GCTA files.
"""

__author__      = "Adriaan van der Graaf"
__copyright__   = "Copyright 2017, Adriaan van der Graaf"


import numpy as np
import scipy
from .. import file_utils
from .. import variants


class CojoCmaFile:
    def __init__(self, file_loc, name):
        self.name = name
        self.ma_results = {}
        with open(file_loc, 'r') as f:
            f.readline()
            for line in f:
                tmp = CojoCmaLine(line)
                self.ma_results[tmp.snp_name] = tmp

    def snps_with_data(self):
        return set(self.ma_results.keys())

    def write_cojo_result(self, file_name):
        lines = [''] * (len(self.ma_results.keys())+1)
        lines[0] = "snp_name\tchr\tbp\tbeta\tse\tp_val\tassoc_name"
        indice = 1
        for i in list(self.ma_results.keys()):
            tmp = self.ma_results[i]
            lines[indice] = '\t'.join([tmp.snp_name,
                               str(tmp.chr),
                               str(tmp.bp),
                               str(tmp.beta_corrected),
                               str(tmp.se_corrected),
                               str(tmp.p_corrected),
                               str(self.name)
                               ]
                              )
            indice += 1


        file_utils.write_list_to_newline_separated_file(lines, file_name)


class CojoCmaLine:
    def __init__(self, line):
        split = [x for x in line[:-1].split() if x != ""]
        self.chr = int(split[0])
        self.snp_name = split[1]
        self.bp = int(split[2])
        self.ref_allele = split[3]
        self.allele_freq = float(split[4])
        self.beta_initial = float(split[5])
        self.se_initial = float(split[6])
        self.p_initial = float(split[7])
        self.freq_geno  = float(split[8])
        self.n_estimated = float(split[9])
        self.beta_corrected = float(split[10])
        self.se_corrected = float(split[11])
        self.se_corrected_original = self.se_corrected
        self.p_corrected = float(split[12])
        self.z_score_corrected = self.beta_corrected / self.se_corrected
        self.has_corrected_se = False

    def get_beta(self):
        return self.beta_corrected

    def get_z_score(self):
        if self.has_corrected_se:
            self.z_score_corrected_ld = self.beta_corrected / self.se_corrected_ld
            return self.z_score_corrected_ld
        else:
            self.z_score_corrected = self.beta_corrected / self.se_corrected
            return self.z_score_corrected

    #todo This method needs to be validated
    def correct_score_based_on_ld(self, r):
        self.has_corrected_se = True
        self.se_corrected_ld = self.se_corrected / abs(r)
        self.z_score_corrected_ld = self.beta_corrected / self.se_corrected_ld



class CojoLdrFile:
    def __init__(self, file_loc, name):
        self.name = name
        with open(file_loc, 'r') as f:
            # first line contains the SNPs
            self.snps = [variants.SNP(x) for x in f.readline()[:-1].split() if (x != '') and (x != 'SNP')]
            self.snp_names = [x.snp_name for x in self.snps]
            self.ld_mat = np.zeros((len(self.snps),len(self.snps)))
            indice = 0
            for line in f:
                self.ld_mat[:,indice] = [float(x) for x in line.split()[1:] if x != '']
                indice += 1

    def get_ld_mat(self):
        return self.ld_mat, self.snp_names

    # this may need to be rewritten, to ensure the snps go from lowest bp position to highest.
    # Not checked right now

    #right now assuming there are no trans effects


    def write_ld_mat_gg(self, filename, bim_data):
        snpnames = [x.snp_name for x in self.snps]

        position = []

        for i in snpnames:
            try:
                position.append(bim_data.bim_results[i].position())
            except:
                position.append(np.NaN)

        ordering = np.argsort(np.array(position))
        string_list = ['chr\tbp\tsnp_name\t' + '\t'.join(np.array(snpnames)[ordering])]
        for i in ordering:
            try:
                tmp = bim_data.bim_results[snpnames[i]].chromosome() + '\t' \
                      + bim_data.bim_results[snpnames[i]].position() + '\t' \
                      + snpnames[i] + '\t'
            except:
                tmp = 'NA\tNA\t' + snpnames[i] + '\t'

            tmp += '\t'.join([str(x) for x in self.ld_mat[i, ordering]])
            string_list.append(tmp)

        file_utils.write_list_to_newline_separated_file(string_list, filename)


class MaFile:
    def __init__(self, file_loc, name):
        self.name = name
        self.ma_results = {}
        with open(file_loc, 'r') as f:
            f.readline()
            for line in f:
                tmp = MaLine(line)
                self.ma_results[tmp.snp_name] = tmp

    def snp_names(self, no_palindromic = False):
        if not no_palindromic:
            return(self.ma_results.keys())
        else:
            palindromic = ["GC", "CG", "AT", "TA"]
            snpnames = []
            for i in self.ma_results.keys():
                if self.ma_results[i].allele_1 + self.ma_results[i].allele_2 not in palindromic:
                    snpnames.append(i)
            return snpnames

    def write_result(self, file_name):
        lines = [''] * (len(self.ma_results.keys())+1)
        lines[0] = "snp_name\tbeta\tse\tp_val\tassoc_name\tbp\tchr"
        indice = 1
        for i in list(self.ma_results.keys()):
            tmp = self.ma_results[i]
            if tmp.has_pos_chr:
                lines[indice] = '\t'.join([tmp.snp_name, str(tmp.beta), str(tmp.se), str(tmp.p_value), self.name,
                                           str(tmp.pos), str(tmp.chr)])
            else:
                lines[indice] = '\t'.join([tmp.snp_name, str(tmp.beta), str(tmp.se), str(tmp.p_value), self.name,
                                          "NA", "NA"])
            indice +=1
        file_utils.write_list_to_newline_separated_file(lines, file_name)

    def add_bim_data(self, bim_data):
        for i in self.ma_results.keys():
            if i in bim_data.bim_results.keys():
                self.ma_results[i].add_pos_chr(bim_data.bim_results[i].position(),
                                               bim_data.bim_results[i].chromosome()
                                               )


class MaLine():
    def __init__(self, line):
        split = [x for x in line[:-1].split() if x != ""]
        self.snp_name = split[0]
        self.allele_1 = split[1]
        self.allele_2 = split[2]
        self.allele_freq = float(split[3])
        self.beta = float(split[4])
        self.se = float(split[5])
        self.p_value = float(split[6])
        self.n_individuals = float(split[7])
        self.has_pos_chr = False

    def add_pos_chr(self, pos, chr):
        self.pos = pos
        self.chr = chr
        self.has_pos_chr = True


def isolate_snps_from_list(snp_loc, gwas_in, gwas_out):
    snp_list = {}

    with open(snp_loc, 'r') as f:
        for line in f:
            snp_list[line[:-1]] = 0

    write_file = open(gwas_out, 'w')

    with open(gwas_in, 'r') as f:
        for line in f:
            split = [x for x in line.split() if x != '']
            if split[0] in snp_list.keys():
                write_file.write(line)

    write_file.close()

    return MaFile(gwas_out, gwas_out)


class SmrMaFile:
    def __init__(self, file_loc, name):

        self.name = name
        self.ma_results = {}
        with open(file_loc, 'r') as f:
            f.readline()
            for line in f:
                tmp = SmrResult(line)
                self.ma_results[tmp.probe_id] = tmp

    def snp_names(self):
        return(self.ma_results.keys())


class SmrResult:
    def __init__(self, line):
        split = [x for x in line[:-1].split() if x != ""]

        self.probe_id = split[0]
        self.probe_chr = int(split[1])
        self.gene = split[2]
        self.probe_bp = int(split[3])
        self.top_snp = split[4]
        self.top_snp_chr = int(split[5])
        self.top_snp_bp = int(split[6])
        self.allele_1 = split[7]
        self.allele_2 = split[8]
        self.allele_freq = float(split[9])
        self.beta_gwas = float(split[10])
        self.se_gwas = float(split[11])
        self.p_gwas = float(split[12])
        self.beta_eqtl = float(split[13])
        self.se_eqtl = float(split[14])
        self.p_eqtl = float(split[15])
        self.beta_smr = float(split[16])
        self.se_smr = float(split[17])
        self.p_smr = float(split[18])
        try:
            self.p_het = float(split[19])
            self.nsnp_het = float(split[20])
            self.has_het_test = True
        except ValueError:
            self.has_het_test = False

    def get_beta_outcome(self):
        beta = self.beta_gwas
        return beta

    def get_z_score_outcome(self):
        return self.beta_gwas / self.se_gwas

    def get_beta_exposure(self):
        return self.beta_eqtl

    def get_z_score_exposure(self):
        return self.beta_gwas / self.se_gwas



def estimate_beta_se_smr(smr_result):
    """
    Determine SMR test effect and standard error.

    :param exposure_data: Exposure estimates which must have the methods get_beta and get_z_score
    :param outcome_data: Outcome estimates which must have the methods get_beta and get_z_score
    :return: tuple of smr beta and smr se
    """

    # from SMR paper.
    t_stat = ((smr_result.get_z_score_exposure() ** 2) * (smr_result.get_z_score_outcome() ** 2)) \
             /  \
             ((smr_result.get_z_score_exposure() ** 2) + (smr_result.get_z_score_outcome() ** 2))

    # this needs a check,
    # checked it with results from Zhy et al. see validate_SMR_estimates.py
    p_value = scipy.stats.chi2.sf(t_stat, 1)
    z_score = scipy.stats.norm.isf(p_value / 2)

    # also from Zhu et al.
    beta_smr = smr_result.get_beta_outcome() / smr_result.get_beta_exposure()
    se_smr = abs(beta_smr / z_score)

    return beta_smr, se_smr, p_value
