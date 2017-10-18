import gzip
import struct
import scipy.stats
from . association_classes import *


def add_p_values_to_associations_dict(associations):
    """
    This class will add wald p values to Association classes.

    :param associations:
    :return:
    """

    assoc_names = list(associations.keys())
    abs_z_scores = [abs(associations[assoc_name].z_score) for assoc_name in assoc_names]
    p_vals = scipy.stats.norm.sf(abs_z_scores) * 2
    [associations[assoc_names[i]].set_wald_p_value(p_vals[i]) for i in range(len(assoc_names))]
    return associations, p_vals


def turn_bin_into_plink_assoc(bytes):
    return struct.unpack('hiifff', bytes)


def read_bin_file(file_name, bim_data):
    # read the gene name:
    with gzip.open(file_name, "rb") as f:
        full_array = f.read()  # should not be so big, maybe a hundred Mb

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
        except RuntimeError:
            continue

        gene_num += 1

    associations, p_values_vec = add_p_values_to_associations_dict(associations)

    return associations, p_values_vec