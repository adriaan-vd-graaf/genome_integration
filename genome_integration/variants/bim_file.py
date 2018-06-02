from genome_integration.variants.SNP import *
import numpy as np

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
