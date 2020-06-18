from genome_integration import utils
from genome_integration import simulate_mr
import numpy as np
import hashlib
from functools import partial

#taken from Stack overflow
def md5sum(filename):
    with open(filename, mode='rb') as f:
        d = hashlib.md5()
        for buf in iter(partial(f.read, 128), b''):
            d.update(buf)
    return d.hexdigest()


def test_plink_file_writing():
    rel_path = '/'.join(('test_resources', 'subset_of_exposure_cohort'))
    if len(__file__.split("/")) > 1:
        plink_loc = "{}/{}".format("/".join(__file__.split("/")[:-1]), rel_path)
    else:
        plink_loc = rel_path

    out_path = '/'.join(('temp_data', 'check_bed_writing'))
    if len(__file__.split("/")) > 1:
        out_loc = "{}/{}".format("/".join(__file__.split("/")[:-1]), out_path)
    else:
        out_loc = out_path


    plinkfile = utils.PlinkFile(plink_loc)
    plinkfile.read_bed_file_into_numpy_array(allele_2_as_zero=False)
    orig_bed_md5sum, orig_bim_md5sum, orig_fam_md5sum = [md5sum(f'{rel_path}{x}') for x in ['.bed', '.bim', '.fam']]

    #this is where all the files are made.
    plinkfile.output_genotypes_to_bed_file(out_loc)

    #fam file comparison
    new_fam_md5sum = md5sum(out_loc + '.fam')
    assert(new_fam_md5sum == orig_fam_md5sum)

    # bim file comparison
    new_bim_md5sum = md5sum(out_loc + '.bim')
    assert(new_bim_md5sum == orig_bim_md5sum)

    #bed_file_comparison
    new_bed_md5sum = md5sum(out_loc + '.bed')

    assert(new_bed_md5sum == orig_bed_md5sum)

if __name__ == '__main__':
    test_plink_file_writing()