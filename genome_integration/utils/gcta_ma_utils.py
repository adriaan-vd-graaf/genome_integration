import numpy as np
from .file_utils import *
from .. import association

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
                if self.ma_results[i].major_allele + self.ma_results[i].minor_allele not in palindromic:
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
        write_list_to_newline_separated_file(lines, file_name)

    def add_bim_data(self, bim_data):
        for i in self.ma_results.keys():
            if i in bim_data.bim_results.keys():
                self.ma_results[i].add_pos_chr(bim_data.bim_results[i].position,
                                               bim_data.bim_results[i].chromosome
                                               )

    #this is ugly.
    def update_snp_names(self):
        tmp_results = {}
        for x in self.ma_results:
            tmp_results[self.ma_results[x].snp_name] = self.ma_results[x]

        self.ma_results = tmp_results

    def delete_everything_except_set(self, snp_set):
        """
        DANGEROUS TO USE, do not use, except if you really want to delete some data in this class

        :param snp_set:
        :returns: itself.
        """
        snp_set = set(snp_set) # make it a set, so I can add lists and stuff.
        temp_dict = {}
        for snp in snp_set:
            if snp in self.ma_results.keys():
                temp_dict[snp] = self.ma_results[snp]

        self.ma_results = temp_dict


class MaLine(association.GeneticAssociation):
    def __init__(self, line):

        split = [x for x in line[:-1].split() if x != ""]
        allele_1 = split[1]
        allele_2 = split[2]
        frq = float(split[3])
        if frq > 0.5:
            major = allele_2
            minor = allele_1
        else:
            major = allele_1
            minor = allele_2

        super().__init__(
            dependent_name="unknown",
            explanatory_name=split[0],
            n_observations=int(split[7]),
            beta=float(split[4]),
            se=float(split[5]),
            r_squared=None,
            chromosome=None,
            position=None,
            major_allele=major,
            minor_allele=minor,
            minor_allele_frequency=frq
        )

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