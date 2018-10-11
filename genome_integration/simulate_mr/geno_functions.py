import numpy as np
from plinkio import plinkfile
import statsmodels.api as sm

def geno_frq(geno_vec):
    return 1 - np.sum(geno_vec) / (2 * geno_vec.shape[0])

def scale_geno_vec(geno_vec):
    freq = np.sum(geno_vec) / (2 * geno_vec.shape[0])
    return geno_vec - 2 * freq

def do_gwas_on_scaled_variants(geno_vec, dependent):
    # geno_vec = sm.add_constant(geno_vec) # this is not what you want.
    geno_vec = scale_geno_vec(geno_vec)
    model = sm.OLS(dependent, geno_vec)
    results = model.fit()
    return results.params[0], results.bse[0]

def read_geno_mat(bed_location):
    plink_file = plinkfile.open(bed_location)

    genotype_mat = np.zeros((len(plink_file.loci),len(plink_file.samples)), dtype=float)
    snp_names = [x.name for x in plink_file.loci]
    i_ids = [x.fid + " " +  x.iid for x in plink_file.samples]

    for i in range(len(plink_file.loci)):
        check_object = plink_file.__next__()
        genotype_mat[i, :] = check_object

    return genotype_mat, plink_file


def harmonize_genotype_matrices(ref_loci, compare_loci, compare_geno):
    """

    Will return a harmonized genotype matrix, based on comparison.

    :param ref_plinkio:
    :param ref_geno:
    :param compare_plinkio:
    :param compare_geno:
    :return: harmonized compare_geno matrix

    """

    # Runtime check, make sure the files are in the same ordering.

    name_differences = np.sum([ref_loci[x].name != compare_loci[x].name for x in range(len(ref_loci))])
    if name_differences != 0:
        raise ValueError("Could not harmonize genotype matrices, as they are composed of different snps")

    allele_differences = np.asarray(["{}/{}".format(ref_loci[x].allele1, ref_loci[x].allele2) ==
                                     "{}/{}".format(compare_loci[x].allele2, compare_loci[x].allele1)
                                     for x in range(len(ref_loci))])

    compare_geno[allele_differences,:] = np.abs(compare_geno[allele_differences,:] - 2)

    return compare_geno

def prune_for_a_region(region, plinkio, geno):
    loci = plinkio.get_loci()
    indices_to_keep = []
    plinkio_loci_to_keep = []
    indice = 0
    for locus in loci:
        if region.snp_in_region(locus.chromosome, locus.bp_position):
            indices_to_keep.append(indice)
            plinkio_loci_to_keep.append(locus)
        indice += 1

    new_geno = geno[np.asarray(indices_to_keep),:]

    return new_geno, plinkio_loci_to_keep


def make_grm_and_ld_from_geno_mat(geno_mat):
    """
    :param geno_mat: N x n numpy matrix, n number of individuals, N number of variants, has values 0, 1, 2
    :return: grm: N x N matrix which is scaled by the number of individuals.
    """

    scaled_geno_mat = np.apply_along_axis(scale_geno_vec, 1, geno_mat)
    grm = None # Not implemented.
    ld = np.corrcoef(geno_mat)

    return grm, ld


if __name__ is '__main__':
    pass