"""
These classes are used to implement some features of variants.
"""

__author__      = "Adriaan van der Graaf"
__copyright__   = "Copyright 2017, Adriaan van der Graaf"


import numpy as np
from .. import file_utils

class SNP:
    def __init__(self,
                 snp_name,
                 chromosome=None,
                 position=None,
                 major_allele =None,
                 minor_allele=None,
                 minor_allele_frequency=None
                 ):

        self.snp_name = snp_name
        self.chromosome = chromosome
        self.position = position
        self.major_allele = major_allele
        self.minor_allele = minor_allele
        self.minor_allele_frequency = minor_allele_frequency

        self.has_position_data = self.position != None and self.chromosome != None
        self.has_allele_data = self.major_allele != None and self.minor_allele != None
        self.has_frequency_data = self.minor_allele_frequency != None

    def add_frequency(self, frq):
        self.minor_allele_frequency = frq
        self.has_frequency_data = True


class BimFile:
    def __init__(self, file_name):
        self.bim_results = {}
        self.bim_results_by_pos = {}
        self.snp_names = []
        self.has_ld_data = False
        with open(file_name, 'r') as f:
            for line in f:
                tmp = BimLine(line)
                posname = str(tmp.snp.chromosome) + ":" + str(tmp.snp.position)
                self.snp_names.append(tmp.snp_name())
                self.bim_results[tmp.snp_name()] = tmp
                self.bim_results_by_pos[posname] = tmp

    def write_bim(self, filename):
        with open(filename, 'w') as f:
            f.write("chr\tsnp_name\tCm\tbp\tAllele_1\tAllele_2\n")
            for i in self.bim_results:
                tmp = self.bim_results[i]
                f.write('\t'.join(tmp.split) + "\n")

    def add_frq_information(self, file_name):
        with open(file_name, 'r') as f:
            f.readline()
            for line in f:
                split = [x for x in line.split() if x != ""]
                snp_name = split[1]
                self.bim_results[snp_name].add_minor_allele_frequency(split[2], split[3], float(split[4]))

    def add_ld_mat(self, file_name):
        self.has_ld_data = True
        self.ld_mat = np.zeros([len(self.snp_names), len(self.snp_names)])
        with open(file_name, 'r') as f:
            i = 0
            for line in f:
                self.ld_mat[i,:] = [float(x) for x in line[:-1].split("\t")]
                i += 1

    #todo perhaps add a method to change the se based on the ld, implemented in the cojoLine object.
    #todo make sure the algorithm is correct, making sure snps in high ld are not chosen twice.
    def isolate_LD_similar_snps(self, snp_list_1, snp_list_2, min_ld):

        if not self.has_ld_data:
            raise RuntimeError("LD data was not found in this bim file.")

        overlapping_snps = list(snp_list_1 & snp_list_2)

        snps_to_check = np.array([i for i in range(len(self.snp_names)) if self.snp_names[i] in snp_list_1])
        snps_to_overlap = np.array([i for i in range(len(self.snp_names)) if self.snp_names[i] in snp_list_2])
        snps_to_overlap = [x for x in snps_to_overlap if self.snp_names[x] not in overlapping_snps]

        best_ld = {}
        # fill it up with the correct snps
        for i in overlapping_snps:
            best_ld[i] = (i, 1.0)

        # now look for cool new snps that we can use for beta SMR.
        for i in snps_to_check:
            if self.snp_names[i] in overlapping_snps:
                continue
            try:
                maximum = max(self.ld_mat[i, snps_to_overlap])
                minimum = min(self.ld_mat[i, snps_to_overlap])
            except:
                continue

            if abs(minimum) > maximum:
                ld_max = minimum
            else:
                ld_max = maximum

            if abs(ld_max) >= min_ld:
                which_max = int(np.where(self.ld_mat[i, snps_to_overlap] == ld_max)[0])
                overlapping_snp = self.snp_names[snps_to_overlap[which_max]]
                overlapping_snps.append(overlapping_snp)  # This may be suboptimal. The best snp may not be chosen, as it may be encountered further down in the loop.
                best_ld[self.snp_names[i]] = (overlapping_snp, ld_max)

        return best_ld

    def write_ld_mat(self, filename):
        snpnames = self.snp_names

        position = []

        for i in snpnames:
            try:
                position.append(self.bim_results[i].position())
            except:
                position.append(np.NaN)

        ordering = np.argsort(np.array(position))
        string_list = ['chr\tbp\tsnp_name\t' + '\t'.join(np.array(snpnames)[ordering])]
        for i in ordering:
            try:
                tmp = self.bim_results[snpnames[i]].chromosome() + '\t' \
                      + self.bim_results[snpnames[i]].position() + '\t' \
                      + snpnames[i] + '\t'
            except:
                tmp = 'NA\tNA\t' + snpnames[i] + '\t'

            tmp += '\t'.join([str(x) for x in self.ld_mat[i, ordering]])
            string_list.append(tmp)

        file_utils.write_list_to_newline_separated_file(string_list, filename)


class BimLine:
    def __init__(self, line):
        self.split = [x for x in line[:-1].split() if x != ""]
        self.line = line
        self.snp = SNP(self.split[1], self.split[0], self.split[3], self.split[4], self.split[5])

    def snp_name(self):
        return self.snp.snp_name

    def position(self):
        return self.snp.position

    def chromosome(self):
        return self.snp.chromosome

    def minor_allele(self):
        return self.snp.minor_allele

    def major_allele(self):
        return self.snp.major_allele

    def minor_allele_frequency(self):
        return self.snp.minor_allele_frequency

    def add_minor_allele_frequency(self,  major, minor, freq):
        if (self.snp.major_allele == major) and (self.snp.minor_allele == minor):
            self.snp.add_frequency(freq)
        elif (self.snp.major_allele == minor) and (self.snp.minor_allele == major):
            self.snp.add_frequency(1 - freq)
        else:
            raise RuntimeError("Alleles do not match in snp" + self.snp_name())