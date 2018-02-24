import scipy.stats
import numpy as np
from .. import variants


class BaseAssociation:
    """
    This class is the base class for any association.
    Currently will only contain the data for a linear association.
    No
    """
    def __init__(self, dependent_name=None,
                 explanatory_name=None,
                 n_observations=None,
                 beta=None,
                 se=None,
                 r_squared=None):

        self.dependent_name = dependent_name
        self.explanatory_name = explanatory_name
        self.beta = beta
        self.se = se
        self.n_observations = n_observations
        self.r_squared = r_squared


class Association(BaseAssociation):
    """
    This class will be any normal linear association.
    """
    __slots__ = ['dependent_name', 'explanatory_name', 'beta', 'se', 'n_observations', 'r_squared', 'z_score',
                 'wald_p_val', 'snp']

    def __init__(self, dependent_name, explanatory_name, n_observations, beta, se, r_squared=None):

        super().__init__(
            dependent_name = dependent_name,
            explanatory_name= explanatory_name,
            beta=float(beta),
            se=float(se),
            n_observations=int(float(n_observations)),
            r_squared=r_squared
        )

        if self.se == 0:
            self.se = np.nextafter(0.0, 1)
            self.z_score = np.sign(beta) * 1337 #big enough. this will introduce some bug in super low p values. but whatever.
        else:
            self.z_score = self.beta / self.se

        self.wald_p_val = None  # not calculating it here, is better if calculation is done later.

        self.snp = None

    def set_wald_p_value(self, pval):
        self.wald_p_val = pval


class GeneticAssociation(Association, variants.BaseSNP):
    """
    This class will represent a genetic association, and will probably be a subclass sometime.

    By definition of this class:

    !!!
    THE MINOR ALLELE IS THE EFFECT ALLELE
    !!!

    A decision Which will probably bite me in the behind when trying to integrate multiallelic snps, but, you know...

    It will contain snp and association date.

    """

    __slots__ = ['snp_name', 'chromosome', 'position', 'major_allele', 'minor_allele', 'minor_allele_frequency',
                 'has_position_data', 'has_allele_data', 'has_frequency_data', 'dependent_name', 'explanatory_name',
                 'beta', 'se', 'n_observations', 'r_squared', 'z_score', 'wald_p_val', 'snp', 'reference_allele',
                 'effect_allele', 'alleles']

    def __init__(self,
                 dependent_name,
                 explanatory_name,
                 n_observations,
                 beta,
                 se,
                 r_squared=None,
                 chromosome=None,
                 position=None,
                 major_allele=None,
                 minor_allele=None,
                 minor_allele_frequency=None,
                 reference_allele=None,
                 effect_allele=None
                 ):

        Association.__init__(self,
                             dependent_name,
                             explanatory_name,
                             n_observations,
                             beta,
                             se,
                             r_squared
                             )

        variants.BaseSNP.__init__(self,
                                  explanatory_name,
                                  chromosome,
                                  position,
                                  major_allele,
                                  minor_allele,
                                  minor_allele_frequency
                                  )

        self.alleles = [self.major_allele, self.minor_allele]

        # ensure the reference alleles are initiated.
        # as well as ensuring that the reference alleles match the major and minor alleles.

        if reference_allele is None:
            self.reference_allele = self.major_allele

        if effect_allele is None:
            self.effect_allele = self.minor_allele



        if (not (reference_allele is None))  and (reference_allele not in self.alleles):
            raise ValueError("Reference allele does not match major or minor allele")

        if (not (effect_allele is None)) and (effect_allele in self.alleles):
            raise ValueError("Effect allele does not match major or minor allele")

    def __str__(self):
        try:
            return "{}-{}, {}/{}, {}, {}, {}".format(self.explanatory_name,
                                                     self.dependent_name,
                                                     self.major_allele,
                                                     self.minor_allele,
                                                     self.beta,
                                                     self.se,
                                                     self.wald_p_val)
        except:
            return "{}-{}, {}/{}, {}, {}".format(self.explanatory_name,
                                                     self.dependent_name,
                                                     self.major_allele,
                                                     self.minor_allele,
                                                     self.beta,
                                                     self.se)





    def add_snp_data(self, snp_data, overwrite=False):
        """
        UNTESTED

        This class will return itself with updated snp data.
        It will only change data from a class if the snp_name is the same, or if the position

        Author comment: This is bloody hard to get right.

        :param snp_data, a baseSNP object or bigger.:
        :return self:
        """

        has_updated_position, has_updated_alleles, alleles_flipped, has_updated_frequency = variants.BaseSNP.add_snp_data(self, snp_data)

        if alleles_flipped:
            self.beta *= -1
            self.z_score *= -1


    def make_gcta_ma_header(self):
        """
        WILL NOT TEST

        Will create an ma header.

        :return: String with an ma file header.
        """
        return "SNP\tA1\tA2\tfreq\tb\tse\tp\tN"


    def make_gcta_ma_line(self):
        """
        WILL NOT TEST

        Makes a GCTA line of the genetic variant.

        Will only return a string, will not write to a file, the user is expected to do this himself.

        :return tab separated string that can be part of ma file:
        """

        # make sure the data that we need is available.
        if not self.has_position_data or not self.has_allele_data or not self.has_frequency_data:
            raise RuntimeError("Cannot write an Ma line. Does not contain the necessary data")

        if self.wald_p_val == None:
            raise RuntimeError("No p value present")

        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.snp_name, self.minor_allele, self.major_allele,
                                                         self.minor_allele_frequency, self.beta, self.se, self.wald_p_val,
                                                         self.n_observations)