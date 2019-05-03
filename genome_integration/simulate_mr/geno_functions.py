import numpy as np
from plinkio import plinkfile
import statsmodels.api as sm

def geno_frq(geno_vec, remove_three_encoding = True, encoding_of_a1=1):
    """
    Return the frequency of a genotype vector

    :param geno_vec: numpy vector of genotypes.
    :param remove_three_encoding: remove three encoding from plinkio
    :return: allele frequency of the allele encoding the one.
    """
    if remove_three_encoding:
        # missing values are encoded as 3
        filter = geno_vec != 3.0
    else:
        filter = np.ones(geno_vec.shape, dtype=bool)
    filtered_vec = geno_vec[filter]
    if encoding_of_a1 == 1:
        return np.sum(filtered_vec) / (2 * filtered_vec.shape[0])
    elif encoding_of_a1 == 0:
        return 1 - np.sum(filtered_vec) / (2 * filtered_vec.shape[0])
    else:
        raise ValueError("Allele 1 encoding should be one (1) or zero (0)")

def scale_geno_vec(geno_vec, remove_three_encoding = True):
    """
    Scale a genotype vector in (0,1,2) to mean zero, variance one.
    based on the allele frequency.

    :param geno_vec: numpy vector of genotypes
    :param remove_three_encoding: remove the three encoding from plinkio, sets values to zero.
    :return: scaled genotype vector.
    """
    if remove_three_encoding:
        # plinkio encodes missing values as 3
        filter = geno_vec != 3.0
    else:
        filter = np.ones(geno_vec.shape, dtype=bool)

    filtered_vec = geno_vec[filter]
    freq = np.sum(filtered_vec) / (2 * filtered_vec.shape[0])
    scaled_geno_vec = geno_vec - 2 * freq
    scaled_geno_vec[~filter] = 0.0
    return scaled_geno_vec


def do_gwas_on_scaled_variants(geno_vec, dependent, remove_three_encoding = True):
    """
    Will do a univariate association on a genotype vector and a phenotype vector.

    :param geno_vec: numpy array of (n): genotypes (0,1,2, (missing:3) )
    :param dependent: numpy array of (n): phenotypes (R)
    :param remove_three_encoding: bool, to remove genotypes that are encoded as three.
    :return: tuple of beta and se estimate.
    """

    geno_vec = scale_geno_vec(geno_vec,
                              remove_three_encoding=remove_three_encoding)
    geno_vec = sm.add_constant(geno_vec)
    model = sm.OLS(dependent, geno_vec)
    results = model.fit()

    return results.params[1], results.bse[1]


def read_geno_mat_plinkio(bed_location):
    """
    Reads in a genotype matrix, using plinkio,
    Takes a lot of mem, and a lot of time if there are a lot of individuals.

    :param bed_location: path of the bed file for reading. make sure you've limited the  it enough before you continue.
    :return: genotype matrix, and the plinkio documentation

    """
    plink_file = plinkfile.open(bed_location)

    genotype_mat = np.zeros((len(plink_file.loci),len(plink_file.samples)), dtype=float)
    snp_names = [x.name for x in plink_file.loci]
    i_ids = [x.fid + " " +  x.iid for x in plink_file.samples]

    for i in range(len(plink_file.loci)):
        check_object = plink_file.__next__()


        genotype_mat[i, :] = check_object

    return genotype_mat, plink_file


def ld_from_geno_mat(geno_mat):
    """
    make an LD mat,
    :param geno_mat: n x m numpy matrix,
        n number of individuals, m number of variants, has values 0, 1, 2 and missing is 3
    :param axis: int
    axis on which to apply the function.
    :return: pearson correlation matrix
    """

    scaled_geno_mat = np.apply_along_axis(scale_geno_vec, 1, geno_mat)
    ld = np.corrcoef(scaled_geno_mat.T)
    return ld


def do_gcta_cojo_conditional(reference_geno, associations, indices_of_interest, ref_loci):
    """

    Recreates the conditional analysis of the GCTA cojo paper.
    Does NOT do the stepwise regression, only does joint effects.

    In testing are accurate up to 3 significant digits due to (I'm assuming) floating point differences.

    See the COJO paper (JIAN YANG 2011, Nature Genetics) supplemental methods for the derivation.
    Their software is also open source, for an implementation, there.


    :param reference_geno:
        Genotype matrix

    :param associations: dict
        Dict of the GeneticAssociation class containing a superset of associations to do cojo from,.

    :param indices_of_interest: nparray of indices.
        Indices of ref_loci with which we want to do associations.

    :param ref_loci: list of str
        The names  of the SNPs which are selected for COJO (keys for associations)

    :return: numpy array of shape (len(indices_of_interest), 2)
        COJO results: conditional beta and se on the columns.


    """

    sample_sizes = np.array([associations[x].n_observations for x in ref_loci], dtype=int)

    beta_se_maf = np.asarray(
                           [[associations[x].beta,
                             associations[x].se,
                             associations[x].minor_allele_frequency
                             ]
                         for x in ref_loci]
                        , dtype=float)

    full_d = 2 * beta_se_maf[:, 2] * (1 - beta_se_maf[:, 2]) * sample_sizes
    s_squared = beta_se_maf[:,1] ** 2

    y_t_y = full_d * s_squared * (sample_sizes - 1) + full_d * beta_se_maf[:,0] **2
    median_y_t_y = np.median(y_t_y)

    adjusted_sample_size = (median_y_t_y / (full_d  * s_squared)) - (beta_se_maf[:,0] **2 / s_squared) + 1

    new_b = np.zeros( (len(indices_of_interest), len(indices_of_interest)), dtype=float)
    new_d = np.zeros((len(indices_of_interest), len(indices_of_interest)), dtype=float)

    w_sums = np.sum(reference_geno[:,indices_of_interest] ** 2, axis=0)

    # now fill in the b and d matrices.
    for j in range(len(indices_of_interest)):
        new_d[j,j] = 2 * beta_se_maf[indices_of_interest[j],2] * \
                    (1 - beta_se_maf[indices_of_interest[j],2]) * \
                    adjusted_sample_size[indices_of_interest[j]]

        #see which sample size is smaller.
        n_vec = np.asarray([x
                            if adjusted_sample_size[indices_of_interest[j]] > x
                            else adjusted_sample_size[indices_of_interest[j]]
                            for x in adjusted_sample_size[indices_of_interest]
                            ], dtype=float)

        p_vec =  2 * beta_se_maf[indices_of_interest[j],2] * \
                    (1 - beta_se_maf[indices_of_interest[j],2]) * \
                    2*beta_se_maf[indices_of_interest,2] * \
                    (1 - beta_se_maf[indices_of_interest,2])

        w_vec = np.sum(reference_geno[:,indices_of_interest[j]] ** 2) * w_sums
        covariances = reference_geno[:,indices_of_interest].transpose() @ reference_geno[:, indices_of_interest[j]]

        new_b[j,:] = n_vec * np.sqrt(  p_vec / w_vec  ) * covariances

    new_b_inv = np.linalg.inv(new_b)
    beta = new_b_inv @ new_d @ beta_se_maf[indices_of_interest, 0]


    return np.asarray([beta, np.sqrt(np.diag(new_b_inv))], dtype=float).transpose()


if __name__ is '__main__':
    pass