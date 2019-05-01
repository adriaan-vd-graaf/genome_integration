import warnings
import scipy.stats
import numpy as np
from ..variants import *
import sys


class BaseAssociation:
    """
    This class will be the summary statistics of a univariate linear association
    without intercept. No data on the association is necessary,
    as this is a parent class of which data should be in a marginal association

    Attributes:
    -----------

    dependent_name: str
        Name of the dependent variable, i.e. the `y` in the solved equation `y=xb+e`

        this name is not used anywhere.

    explanatory_name: str
        Name of the explanatory variable, i.e. the `x` in the solved equation `y=xb+e`

    beta: float or castable to float
        value of the slope variable, i.e. the `b` in the solved equation `y=xb+e`

    se: float or castable to float
        value of the standard error of the slope variable, i.e. the `se(b)` in the solved equation `y=xb+e`

    n_observations: int or castable to int
        The number of observations upon which the equation `y=xb+e` is solved.

    r_squared:  float or castable to float
        The coefficient of determination or how much variance the model explains out of total variance.

    Methods
    -------

    None

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
    This class will be the summary statistics of a univariate linear association
    without intercept. Data on the association is necessary, in contrast to the parent class BaseAssociation

    Attributes:
    -----------

    dependent_name: str
        Name of the dependent variable, i.e. the `y` in the solved equation `y=xb+e`

        this name is not used anywhere.

    explanatory_name: str
        Name of the explanatory variable, i.e. the `x` in the solved equation `y=xb+e`

    beta: float or castable to float
        value of the slope variable, i.e. the `b` in the solved equation `y=xb+e`

    se: float or castable to float
        value of the standard error of the slope variable, i.e. the `se(b)` in the solved equation `y=xb+e`

    n_observations: int or castable to int
        The number of observations upon which the equation `y=xb+e` is solved.

    r_squared:  float or castable to float
        The coefficient of determination or how much variance the model explains out of total variance.

    Methods
    -------

    set_wald_p_val():
        Sets the wald test p value of the estimate.
        Warning:
        This identifying this p value is only sufficient if you have sufficient observations.
        otherwise a t statistic is more meaningful.




    """

    def __init__(self, dependent_name, explanatory_name, n_observations, beta, se, r_squared=None):
        """
        Inits the Association class.

        :param dependent_name:
        :param explanatory_name:
        :param n_observations:
        :param beta:
        :param se:
        :param r_squared:
        """
        super().__init__(
            dependent_name = dependent_name,
            explanatory_name= explanatory_name,
            beta=float(beta),
            se=float(se),
            n_observations=int(float(n_observations)),
            r_squared=r_squared
        )

        if self.se == 0:
            warnings.filterwarnings("ignore")
            self.se = np.nextafter(0.0, 1)
            self.z_score = self.beta / self.se
            warnings.filterwarnings("default")
            print(f"Genetic association '{self.dependent_name}, {self.explanatory_name}': "
                  f"se was zero, z score will be infinite.", file=sys.stderr)

        else:
            self.z_score = self.beta / self.se

        self.wald_p_val = None  # not calculating it here, is better if calculation is done later.
        self.p_val = None

        self.snp = None

    def set_p_val(self, pval):
        """
        This class will set a p value you calculated yourself.

        :param pval:
        :return: self
        """
        
        self.wald_p_val = pval
        self.p_val = pval

    def set_wald_p_val(self, pval):
        """
        This class will set a p value you calculated yourself.

        :param pval:
        :return: self
        """

        warnings.warn("this method has been renamed `set_p_val`, please use this, as it is not strictly a wald p value",
                      DeprecationWarning)

        if not np.isfinite(self.z_score):
            self.wald_p_val = 0.0
        else:
            self.wald_p_val = pval


class GeneticAssociation(Association, SNP):
    """
    This class will represent a genetic association. Depends on the SNP and association parent classes.

    By definition of this class:

    definition: THE MINOR ALLELE IS THE EFFECT ALLELE
    This decicion is not great, but it's grandfathered in.
    but for now it's good enough


    Attributes:
    -----------

    Attributes specific to the GeneticAssociation class.

    effect_allele: str
        Which allele is used as the effect allele, meaning the allele which increases the variable `x`
        Important to know, as then it's possible to identify the risk allele.

    Other attributes are inherited from their respective classes, so they are subject to change.


    Methods
    -------

    Specific to the GeneticAssociation class:

    add_snp_data(self, snp_data, overwrite=False)
        Adds data from another SNP class-like object.
        Masked from the SNP class, as it will also the '-1*b' if alleles are flipped.
        Otherwise will do the same.


    Other methods are inherited from the other classes. See their documentation, as it is subject to change.

    """

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

        SNP.__init__(self,
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

        if type(self.wald_p_val) is float:
            return "{}-{}, {}/{}, {}, {}, {}".format(self.explanatory_name,
                                                     self.dependent_name,
                                                     self.major_allele,
                                                     self.minor_allele,
                                                     self.beta,
                                                     self.se,
                                                     self.wald_p_val)
        else:
            return "{}-{}, {}/{}, {}, {}".format(self.explanatory_name,
                                                     self.dependent_name,
                                                     self.major_allele,
                                                     self.minor_allele,
                                                     self.beta,
                                                     self.se)





    def add_snp_data(self, snp_data, overwrite=False):
        """

        This class will return itself with updated snp data.
        It will only change data from a class if the snp_name is the same, or if the position


        :param snp_data, a SNP object or something that was extended from it:
        :return self: but with updated allele, name information.
        """

        has_updated_position, has_updated_alleles, alleles_flipped, has_updated_frequency \
            = SNP.add_snp_data(self,
                               snp_data,
                               overwrite=overwrite)

        if alleles_flipped:
            self.beta *= -1
            self.z_score *= -1

