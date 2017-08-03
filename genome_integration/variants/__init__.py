"""
These classes are used to implement some features of variants.
"""

__author__      = "Adriaan van der Graaf"
__copyright__   = "Copyright 2017, Adriaan van der Graaf"


import numpy as np
from .. import file_utils

class BaseSNP:
    """

    This is the base SNP class.
    Inherit from this class. If you want to use this class,
    please use SNP, as it is just a wrapper.

    """

    def __init__(self,
                 snp_name,
                 chromosome=None,
                 position=None,
                 major_allele=None,
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


    def add_snp_data(self, snp_data, overwrite=False):

        """
        UNTESTED

        This class will return itself with updated snp data.
        It will only change data from a class if the snp_name is the same, or if the position

        Author comment: This is bloody hard to get right.

        :param snp_data, a baseSNP object or bigger.:
        :return self:

        """

        # runtime checks.
        if snp_data.snp_name != self.snp_name or snp_data.snp_name != (
                        str(snp_data.chromosome) + ":" + str(snp_data.position)):
            raise RuntimeError(
                "No match in SNPs between Association: " + self.snp_name + " and " + snp_data.snp_name)

        if overwrite:
            print("Overwriting snp data, effect directions may be lost.")

        #assuming you want the name ot be changed if you do this.
        if not snp_data.snp_name != self.snp_name:
            self.snp_name == snp_data.snp_name

        if (not self.has_position_data) or overwrite:
            self.position = snp_data.position
            self.chromosome = snp_data.position
            # update the presence of position_data
            self.has_position_data = self.position != None and self.chromosome != None

        elif str(self.chromosome) + ":" + str(self.position) != \
                                str(snp_data.chromosome) + ":" + str(snp_data.position):
            # if there is already data, make sure that the positions are the same.
            raise RuntimeError(
                "No position match in SNPs between Association: " + self.snp_name + " and " + snp_data.snp_name)

        swapped = False

        # get the alleles right, takes more logic that I really wanted.
        if (not self.has_allele_data) or overwrite:

            # make sure there is no funny allele swaps if there is information for one allele.
            # if there are allele swaps. then swap the alleles in the data that is passed to the function.
            if (self.major_allele != None or self.minor_allele != None) and not overwrite:
                # there is information for the major alle   le.
                raise RuntimeError("A SNP with a single allele present is being updated, has not been implemented")

            # elif self.minor_allele != None and not overwrite:
            #     # there is information for the minor allele
            #     if self.minor_allele != snp_data.minor_allele:
            #         # make the switch of alleles
            #         tmp = snp_data.major_allele
            #         snp_data.major_allele = snp_data.minor_allele
            #         snp_data.minor_allele = tmp
            #         swapped = True

            # save it up.
            self.major_allele = snp_data.major_allele
            self.minor_allele = snp_data.minor_allele
            self.has_allele_data = self.major_allele != None and self.minor_allele != None

        elif self.major_allele == snp_data.minor_allele and self.minor_allele == snp_data.major_allele_allele:
            # if there is an allele swap, change the swapped to true, so that the data is there.
            swapped = True

        elif self.major_allele != snp_data.major_allele or self.minor_allele != snp_data.minor_allele:
            raise RuntimeError(
                "No allele match in SNPs between Association: " + self.snp_name + " and " + snp_data.snp_name)

        # Because the last checks made sure the alleles are right (let's hope) I can just change the alleles.
        if (not self.has_frequency_data) or overwrite:
            if swapped:
                snp_data.minor_allele_frequency -= 1

            self.minor_allele_frequency = snp_data.minor_allele_frequency



    def add_pos_chr(self, pos, chr):
        self.position = pos
        self.chromosome = chr
        self.has_position_data = True


    def add_minor_allele_frequency(self,  major, minor, freq):
        #if there are no alleles, then just use this.
        if not self.has_allele_data:
            self.minor_allele_frequency = freq
            self.has_allele_data = True
            return

        #allele data is present, so we need to check what it is.
        if (self.major_allele == major) and (self.minor_allele == minor):
            self.minor_allele_frequency = freq
            self.has_frequency_data = True
        elif (self.major_allele == minor) and (self.minor_allele == major):
            self.minor_allele_frequency = 1 - freq
            self.has_frequency_data = True
        else:
            raise RuntimeError("Alleles do not match in snp" + self.snp_name)

        return

class SNP(BaseSNP):
    #I use this, so that I can use the __slots__ function here.
    __slots__ = ['snp_name', 'chromosome', 'position', 'major_allele', 'minor_allele', 'minor_allele_frequency']
    def __init__(self,
                 snp_name,
                 chromosome=None,
                 position=None,
                 major_allele=None,
                 minor_allele=None,
                 minor_allele_frequency=None
                 ):
        super().__init__(snp_name, chromosome, position, major_allele, minor_allele, minor_allele_frequency)


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

    # Todo, This may be ripe for removal as I'll ensure this method not necessary anymore
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

    # Todo, This may be ripe for removal as I'll ensure this method not necessary anymore
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
        super().__init__(split[1], split[0], split[3], split[4], split[5])
