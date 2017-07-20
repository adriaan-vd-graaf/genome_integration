import numpy as np
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

