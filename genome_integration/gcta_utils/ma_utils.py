
from .. import file_utils

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

    def delete_everything_except_set(self, snp_set):
        """
        DANGEROUS TO USE, do not use, except if you really want to delete some data.

        :param snp_set:
        :return:
        """
        snp_set = set(snp_set) # make it a set, so I can add lists and stuff.
        temp_dict = {}
        for snp in snp_set:
            if snp in self.ma_results.keys():
                temp_dict[snp] = self.ma_results[snp]

        self.ma_results = temp_dict


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

    def get_beta(self):
        return self.beta

    def get_se(self):
        return self.se

    def get_z_score(self):
        return self.beta / self.se

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
