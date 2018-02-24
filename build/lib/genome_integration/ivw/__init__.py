"""
This class is used to do a 'simple'
inverse variance estimate on independent effects.
"""

__author__      = "Adriaan van der Graaf"
__copyright__   = "Copyright 2017, Adriaan van der Graaf"

from . IVWResult import *

# todo, make this more compatible with the IVWResult class.
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
        self.gene_name = split[1]
        self.chromosome = split[2]
        self.n_snps = int(split[3])
        self.pearson_r = float(split[4])
        self.pearson_r_p = float(split[5])
        self.beta_ivw = float(split[6])
        self.se_ivw = float(split[7])
        self.p_ivw = float(split[8])
        self.snp_names = set(split[9].split(","))
        self.snp_smr_beta = split[10].split(",")
        self.snp_smr_se = split[11].split(",")