import scipy
from .. import variants



class Association:
    __slots__ = ['dependent_name', 'explanatory_name', 'beta', 'se', 'n_observations', 'r_squared', 'z_score',
                 'wald_p_val', 'snp']

    def __init__(self, dependent_name, explanatory_name, n_observations, beta, se, r_squared=None):
        self.dependent_name = dependent_name
        self.explanatory_name = explanatory_name
        self.beta = float(beta)
        self.se = float(se)
        self.n_observations = int(float(n_observations))
        self.r_squared = r_squared

        self.z_score = self.beta / self.se #doing this ensures that beta and se are not none.

        self.wald_p_val = None  # not calculating it here, is better if calculation is done later.

        self.snp = None


    def set_wald_p_value(self, p_val):  # this may seem unintuitive, but It's faster to do this outside an array.
        self.wald_p_val = p_val

    def get_beta(self):
        return self.beta

    def get_se(self):
        return self.se

    def get_z_score(self):
        return self.z_score


class GeneticAssociation(Association, variants.BaseSNP):
    """
    This class will represent a genetic association, and will probably be a subclass sometime.

    By definition of this class:

    !!!
    the MINOR ALLELE IS THE EFFECT ALLELE
    !!!

    A decision Which will probably bite me in the behind when trying to integrate multiallelic snps, but, you know...

    It will contain snp and association date.

    """

    __slots__ = ['snp_name', 'chromosome', 'position', 'major_allele', 'minor_allele', 'minor_allele_frequency',
                 'has_position_data', 'has_allele_data', 'has_frequency_data', 'dependent_name', 'explanatory_name',
                 'beta', 'se', 'n_observations', 'r_squared', 'z_score', 'wald_p_val', 'snp']

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
                 minor_allele_frequency=None
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

        if (not self.snp_name == snp_data.snp_name)  and \
                (not self.snp_name == str(snp_data.chromosome) + ":" + str(snp_data.position)):

            raise RuntimeError(
                "No match in SNPs between Association: " + self.snp_name + " and " +
                str(snp_data.chromosome) + ":" + str(snp_data.position))

        if overwrite:
            print("Overwriting snp data, effect directions may be lost.")

        #if the snp_name is a position, update it to an rs number.
        if self.snp_name != snp_data.snp_name:
            self.snp_name = snp_data.snp_name

        if (not self.has_position_data) or overwrite:
            self.position = snp_data.position
            self.chromosome = snp_data.chromosome
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

        elif self.major_allele == snp_data.minor_allele and self.minor_allele == snp_data.major_allele:
            # if there is an allele swap, change the swapped to true, so that the data is there.
            swapped = True

        elif self.major_allele != snp_data.major_allele or self.minor_allele != snp_data.minor_allele:
            raise RuntimeError(
                "No allele match in SNPs between Association: " + self.snp_name +
                " {}/{} and ".format(self.major_allele, self.minor_allele) + snp_data.snp_name +
                " {}/{}.".format(snp_data.major_allele, snp_data.minor_allele) )

        # Because the last checks made sure the alleles are right I can just change the alleles.
        if (not self.has_frequency_data) or overwrite:
            if swapped:
                snp_data.minor_allele_frequency *= -1
                self.beta *= -1
                self.z_score *= -1

            self.minor_allele_frequency = snp_data.minor_allele_frequency
            self.has_frequency_data = True

    def make_gcta_ma_header(self):
        return "SNP\tA1\tA2\tfreq\tb\tse\tp\tN"

    def make_gcta_ma_line(self):
        """
        Makes a GCTA line of the genetic variant.

        Will only submit a string, will not write to a file, the user is expected to do this himself.

        I'm assuming that there are

        :return tab separated string that can be part of ma file:
        """

        # make sure the data that we need is available.
        if not self.has_position_data or not self.has_allele_data or not self.has_frequency_data:
            raise RuntimeError("Cannot write an Ma line. Does not contain the necessary data")

        if self.wald_p_val == None:
            raise RuntimeError("No p value present")

        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.snp_name, self.major_allele, self.minor_allele,
                                                         self.minor_allele_frequency, self.beta, self.se, self.wald_p_val,
                                                         self.n_observations)


    def do_smr_estimate_using_self_as_exposure(self, outcome_Genetic_association):
        """
            Determine SMR test effect and standard error.
        """

        z_score_exposure = self.z_score
        z_score_outcome = outcome_Genetic_association.z_score


        # from SMR paper.
        t_stat = ((z_score_exposure ** 2) * (z_score_outcome ** 2)) \
                 / \
                 ((z_score_exposure ** 2) + (z_score_outcome ** 2))

        # this needs a check,
        # checked it with results from Zhy et al. see validate_SMR_estimates.py
        p_value = scipy.stats.chi2.sf(t_stat, 1)
        z_score = scipy.stats.norm.isf(p_value / 2)

        # also from Zhu et al.
        beta_smr = outcome_Genetic_association.beta / self.beta
        se_smr = abs(beta_smr / z_score)

        return beta_smr, se_smr, p_value

