"""
This class is used to read in eqtl mapping pipeline formats.
"""

import numpy as np
import gzip
from .. association import *

class EQTLMappingAssociation(GeneticAssociation):
    """
    Read in an eQTL mapping pipeline line, and pour it into the geneticassociation class.
    """

    def __init__(self, line):
        split = line.split()

        minor = split[9]
        alleles = split[8].split("/")

        if alleles[0] == minor:
            major = alleles[1]
        else:
            major = alleles[0]

        super().__init__(
            dependent_name=split[16],
            explanatory_name=split[1],
            n_observations=3800,
            beta=np.nan,
            se=np.nan,
            r_squared=np.nan,
            chromosome=split[2],
            position=split[3],
            major_allele=major,
            minor_allele=minor,
            minor_allele_frequency=None
        )

        self.wald_p_val = float(split[0])


def ReadEQTLMappingFile(filelocation):

    gene_dict = {}

    with gzip.open(filelocation, "rb") as f:
        f.readline()
        for line in f:
            line = line.decode("utf8")

            association = EQTLMappingAssociation(line)
            assoc_dict = {association.snp_name: association}
            try:
                gene_dict[association.dependent_name].update(assoc_dict)
            except KeyError:
                gene_dict[association.dependent_name] = assoc_dict

    return gene_dict

