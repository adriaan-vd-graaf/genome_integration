from genome_integration import utils
from genome_integration import simulate_mr
import numpy as np


def test_plink_bed_file_reading():
    """
    Compares this implementation of plink file reading to that of plinkio.

    :return: None, would raise assertionerror if the data is not exactly the same.
    """
    rel_path = '/'.join(('test_resources', 'subset_of_exposure_cohort'))

    if len(__file__.split("/")) > 1:
        plink_loc = "{}/{}".format("/".join(__file__.split("/")[:-1]), rel_path)
    else:
        plink_loc = rel_path

    plinkio = simulate_mr.read_geno_mat_plinkio(plink_loc)

    plinkfile = utils.PlinkFile(plink_loc)
    own_implementation = plinkfile.read_bed_file_into_numpy_array(allele_2_as_zero=False)

    assert(np.all(plinkio[0].T == own_implementation))

    plinkio_snp_names = [x.name for x in plinkio[1].loci]
    plinkio_sample_names = [f"{x.fid}~__~{x.iid}"for x in plinkio[1].samples]

    these_snp_names = [plinkfile.bim_data.bim_results[x].snp_name for x in plinkfile.bim_data.snp_names]
    these_fam_names = [x for x in plinkfile.fam_data.sample_names]

    assert(plinkio_snp_names == these_snp_names)
    assert(plinkio_sample_names == these_fam_names)

    #make sure flipped encoding is correct.
    flipped_encoding = plinkfile.read_bed_file_into_numpy_array(allele_2_as_zero=True)
    plinkio_mat = plinkio[0].T
    assert(np.all((plinkio_mat == 2) == (flipped_encoding == 0)))
    assert(np.all((plinkio_mat == 0) == (flipped_encoding == 2)))
    assert (np.all((plinkio_mat == 3) == (flipped_encoding == 3)))
    assert (np.all((plinkio_mat == 1) == (flipped_encoding == 1)))





test_plink_bed_file_reading()

