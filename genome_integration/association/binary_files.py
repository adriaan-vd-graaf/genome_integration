import gzip
import struct
import scipy.stats
from . association_classes import *

"""
These classes are helper functions to [semi] efficiently store and read association statistics.

The format of the association file and struct package encoding reprsents the following.

The first bytes, delimited by the b'\0' byte are the name of the gene

Then all associations are represented in the following manner until the file ends.
    short (h): chromosome (only numeric chromosome allowed)  
    int (i): position 
    int (i): number of observations
    float (f): `b` slope estimate of the equation `y = xb + e`.
                `y` is phenotype vector, `x` is genotype vector of a variant, `e` is residual error
    float (f): se(b) standard error of the beta estimate
    float (f): R^2 variance explained by the model.
    
    Effect alleles are not stored. These need to be known.
"""

def turn_plink_assoc_line_into_bin(line):
    """
    This turns a non header line from a standard plink repeated univariate
    association file (*.qassoc) into a byte format.
    The format represents in order:

        short (h): chromosome (only numeric chromosome allowed)
        int (i): position
        int (i): number of observations
        float (f): `b` slope estimate of the equation `y = xb + e`.
                    `y` is phenotype vector, `x` is genotype vector of a variant, `e` is residual error
        float (f): se(b) standard error of the beta estimate
        float (f): R^2 variance explained by the model.

    :param line:
    :return: either a byte array with the format, or if something went wrong, an empty byte array.
    """
    try:
        split = [x for x in line.split() if x != ""]
        bytes = struct.pack('hiifff',int(split[0]), int(split[2]), int(split[3]), float(split[4]), float(split[5]), float(split[6]))
        return bytes
    except:
        return b''


def turn_plink_assoc_into_bin(in_file, out_file, gene_name):
    """
    This function parses a plink association file and turns it into a binary association file.

    :param in_file: File name of the plink .qassoc file
    :param out_file: File name to output to.
    :param gene_name: name of the gene for which the associations are made.
    :return:
    """

    outfile = open(out_file, "wb")
    outfile.write(bytearray(gene_name, 'utf8') + b'\0')
    with open(in_file, 'r') as f:
        f.readline() # header
        for line in f:
            outfile.write(turn_plink_assoc_line_into_bin(line))
    outfile.close()




def add_p_values_to_associations_dict(associations):

    """
    This class will add wald p values to Association classes.
    As this is expensive to estimate one at a time, this is done in vector form.

    :param associations:
    :return:
    """

    assoc_names = list(associations.keys())
    abs_z_scores = [abs(associations[assoc_name].z_score) for assoc_name in assoc_names]
    n_obs = np.asarray([abs(associations[assoc_name].n_observations) for assoc_name in assoc_names], dtype=int)
    p_vals = scipy.stats.t.sf(abs_z_scores, n_obs - 2) * 2
    [associations[assoc_names[i]].set_p_val(p_vals[i]) for i in range(len(assoc_names))]
    return associations, p_vals


def turn_bin_into_plink_assoc(bytes):
    """
    Turns a byte array into a tuple of associations.

    :param bytes:
    :return: tuple of ints and floats, representing the summary statistics.
    """
    return struct.unpack('hiifff', bytes)


def read_bin_file(file_name, bim_data):
    """
    Turns a compressed file of associations into a dict that fully contains the GeneticAssociation class.


    :param file_name: File where the binary associations are stored
    :param bim_data: the BimFile object (see genome_integration.variants) that is
                     used to know which alleles were used for effect alele.
    :return: dictionary of all Genetic associations
    """
    # read the gene name:
    with gzip.open(file_name, "rb") as f:
        full_array = f.read()  # should not be so big, maybe a hundred Mb for about 8 mil SNPs

    gene_name, sep, eqtl_data = full_array.partition(b'\0')  # sep is the separator.

    gene_num = 0

    eqtl_size = 24

    associations = {}

    for i in range(int(len(eqtl_data) / eqtl_size)):
        data = turn_bin_into_plink_assoc(eqtl_data[i * eqtl_size:(i * eqtl_size) + eqtl_size])
        pos_name = str(data[0]) + ":" + str(data[1])
        association = GeneticAssociation(gene_name, pos_name, data[2], data[3], data[4], data[5], chromosome=data[0], position=data[1])

        try:
            association.add_snp_data(bim_data.bim_results_by_pos[pos_name])
            association.snp_name = bim_data.bim_results_by_pos[pos_name].snp_name
            associations[pos_name] = association
        except KeyError:  # No SNP data present for the variant, do nothing.
            continue
        except RuntimeError: # No correct data present for the variant, do nothing.
            continue

        gene_num += 1

    associations, _ = add_p_values_to_associations_dict(associations)

    return associations


def read_bin_file_region(file_name, bim_data, gene_region):
    """
    Turns a compressed file of associations into a dict that fully contains the GeneticAssociation class.

    :param file_name: File where the binary associations are stored
    :param bim_data: the BimFile object (see genome_integration.variants) that is
                     used to know which alleles were used for effect alele.
    :param gene_region: The StartEndRegion class, which is used to only filter for a specific region,
                        making it a bit faster and less expensive in mem to read than the full file.

    :return: dictionary of all Genetic associations
    """
    # read the gene name:
    with gzip.open(file_name, "rb") as f:
        full_array = f.read()  # should not be so big, maybe a hundred Mb

    gene_name, sep, eqtl_data = full_array.partition(b'\0')  # sep is the separator.

    gene_num = 0

    eqtl_size = 24

    associations = {}

    for i in range(int(len(eqtl_data) / eqtl_size)):
        data = turn_bin_into_plink_assoc(eqtl_data[i * eqtl_size:(i * eqtl_size) + eqtl_size])

        if not gene_region.snp_in_region(data[0], data[1]):
            continue

        pos_name = str(data[0]) + ":" + str(data[1])
        association = GeneticAssociation(gene_name, pos_name, data[2], data[3], data[4], data[5], chromosome=data[0], position=data[1])

        try:
            association.add_snp_data(bim_data.bim_results_by_pos[pos_name])
            association.snp_name = bim_data.bim_results_by_pos[pos_name].snp_name
            associations[pos_name] = association
        except KeyError:  # No SNP data present for the variant, do nothing.
            continue
        except RuntimeError: # No correct data present for the variant, do nothing.
            continue

        gene_num += 1

    associations, p_values_vec = add_p_values_to_associations_dict(associations)

    return associations, p_values_vec