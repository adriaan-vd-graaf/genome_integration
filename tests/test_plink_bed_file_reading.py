from genome_integration import utils
from genome_integration import simulate_mr
import numpy as np
import time


def test_plink_bed_file_reading():
    """
    Compares this implementation of plink file reading to that of plinkio.

    :return:
    """
    rel_path = '/'.join(('test_resources', 'subset_of_exposure_cohort'))

    plink_loc = "{}/{}".format("/".join(__file__.split("/")[:-1]), rel_path)


    start = time.time()
    plinkio = simulate_mr.read_geno_mat(plink_loc)
    end = time.time()

    start = time.time()
    plinkfile = utils.PlinkFile(plink_loc)
    own_implementation = plinkfile.read_bed_file_into_numpy_array()
    end = time.time()

    check = np.asarray([own_implementation[:,0], plinkio[0][0,:]])
    assert(np.all(plinkio[0].T == own_implementation))

    plinkio_snp_names = [x.name for x in plinkio[1].loci]
    plinkio_sample_names = [f"{x.fid}~__~{x.iid}"for x in plinkio[1].samples]

    these_snp_names = [plinkfile.bim_data.bim_results[x].snp_name for x in plinkfile.bim_data.snp_names]
    these_fam_names = [x for x in plinkfile.fam_data.sample_names]

    assert(plinkio_snp_names == these_snp_names)
    assert(plinkio_sample_names == these_fam_names)



test_plink_bed_file_reading()

