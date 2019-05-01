import numpy as np
from .file_utils import *
from .. import association

class MaFile:
    """
    Implements a reader and writer to the MA file as described by GCTA COJO

    Attributes
    ----------
    name: str
        name of the file.
    ma_results: dict
        keys are the SNP names and values are of class MaLine

    Methods
    -------
    snp_names(self, no_palindromic = False)
        return the snp names in the mafile

    write_results(self, file_name)
        write the results to a file name

    add_bim_data(self, bim_data)
        add the position and chromosome from a BimFile reference to all the SNPS.

    """
    def __init__(self, file_loc, name):
        self.name = name
        self.ma_results = {}
        with open(file_loc, 'r') as f:
            f.readline()
            for line in f:
                tmp = MaLine(line)
                self.ma_results[tmp.snp_name] = tmp

    def snp_names(self, no_palindromic = False):
        """
                return the snp names in the mafile

        :param no_palindromic: bool
            If also palindromic SNPs should be returned (default is False)
        :return:
            list of SNP names.
        """
        if not no_palindromic:
            return(list(self.ma_results.keys()))
        else:
            palindromic = ["GC", "CG", "AT", "TA"]
            snpnames = []
            for i in self.ma_results.keys():
                if self.ma_results[i].major_allele + self.ma_results[i].minor_allele not in palindromic:
                    snpnames.append(i)
            return snpnames

    def write_result(self, file_name):
        """

        write the results to a file name
        :param file_name: str
            file name to write to
        :return: None

        """
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
        """
        add the position and chromosome from a BimFile reference to all the SNPS.

        :param bim_data:
        :return: None
        """
        for i in self.ma_results.keys():
            if i in bim_data.bim_results.keys():
                self.ma_results[i].add_pos_chr(bim_data.bim_results[i].position,
                                               bim_data.bim_results[i].chromosome
                                               )




class MaLine(association.GeneticAssociation):
    """
    Class that inherits from GeneticAssociation, reading an MA line.

    Attributes
    ----------
    All the attributes of the GeneticAssociation class.

    """
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
