"""
These classes are used to implement some features of variants.
"""

__author__      = "Adriaan van der Graaf"
__copyright__   = "Copyright 2017, Adriaan van der Graaf"


import numpy as np

from .. SNP import *
from .. import file_utils


class BimFile:
    def __init__(self, file_name):
        self.bim_results = {}
        self.bim_results_by_pos = {}
        self.snp_names = []
        self.has_ld_data = False
        with open(file_name, 'r') as f:
            for line in f:
                tmp = BimLine(line)
                posname = str(tmp.chromosome) + ":" + str(tmp.position)
                self.snp_names.append(tmp.snp_name)
                self.bim_results[tmp.snp_name] = tmp
                self.bim_results_by_pos[posname] = tmp

    def write_bim(self, filename):
        with open(filename, 'w') as f:
            f.write("chr\tsnp_name\tCm\tbp\tAllele_1\tAllele_2\n")
            for i in self.bim_results:
                tmp = self.bim_results[i]
                f.write('\t'.join(tmp.split) + "\n")


    def add_frq_information(self, file_name):
        """

        :param file_name: file name of the plink frq or frqx file.
        :return: self, with added frq information

        """
        with open(file_name, 'r') as f:
            split = f.readline()[:-1].split("\t")
            # frq file.
            if len(split) == 5:
                for line in f:
                    split = [x for x in line.split() if x != ""]
                    snp_name = split[1]
                    try:
                        self.bim_results[snp_name].add_minor_allele_frequency(split[3], split[2], float(split[4]))
                    except KeyError:
                        try:
                            self.bim_results_by_pos[snp_name].add_minor_allele_frequency(split[3], split[2], float(split[4]))
                        except KeyError:
                            continue
            # frqx file
            elif len(split) == 10:
                for line in f:
                    split = [x for x in line.split() if x != ""]
                    snp_name = split[1]
                    a_one_count = int(split[4])*2 + int(split[5])
                    a_two_count = int(split[6])*2 + int(split[5])
                    if a_one_count <= a_two_count:
                        minor = split[2]
                        major = split[3]
                        maf = float(a_one_count) / float(a_one_count + a_two_count)

                    else:
                        minor = split[3]
                        major = split[2]

                        maf = float(a_two_count) / float((a_one_count + a_two_count))
                    try:
                        self.bim_results[snp_name].add_minor_allele_frequency(major, minor, float(maf))
                    except KeyError:
                        try:
                            self.bim_results_by_pos[snp_name].add_minor_allele_frequency(major, minor, float(maf))
                        except:
                            continue
            else:
                RuntimeError("The frq file header was not in any correct formatting.")


    def add_ld_mat(self, file_name):
        self.has_ld_data = True
        self.ld_mat = np.zeros([len(self.snp_names), len(self.snp_names)])
        with open(file_name, 'r') as f:
            i = 0
            for line in f:
                self.ld_mat[i,:] = [float(x) for x in line[:-1].split("\t")]
                i += 1

    #todo, This may be ripe for removal as I'll ensure this method not necessary anymore
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

    # todo, This may be ripe for removal as I'll ensure this method not necessary anymore
    def write_ld_mat(self, filename):
        snpnames = self.snp_names

        position = []

        for i in snpnames:
            try:
                position.append(self.bim_results[i].position)
            except:
                position.append(np.NaN)

        ordering = np.argsort(np.array(position))
        string_list = ['chr\tbp\tsnp_name\t' + '\t'.join(np.array(snpnames)[ordering])]
        for i in ordering:
            try:
                tmp = self.bim_results[snpnames[i]].chromosome + '\t' \
                      + self.bim_results[snpnames[i]].position + '\t' \
                      + snpnames[i] + '\t'
            except:
                tmp = 'NA\tNA\t' + snpnames[i] + '\t'

            tmp += '\t'.join([str(x) for x in self.ld_mat[i, ordering]])
            string_list.append(tmp)

        file_utils.write_list_to_newline_separated_file(string_list, filename)


class BimLine(BaseSNP):
    __slots__ = ['snp_name', 'chromosome', 'position', 'major_allele', 'minor_allele', 'minor_allele_frequency']

    def __init__(self, line):
        split = [x for x in line[:-1].split() if x != ""]
        super().__init__(snp_name=split[1],
                 chromosome=split[0],
                 position=split[3],
                 major_allele=split[5],
                 minor_allele=split[4],
                 minor_allele_frequency=None
                )
