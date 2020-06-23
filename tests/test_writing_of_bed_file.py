from genome_integration import utils
from genome_integration import gene_regions

import subprocess
import hashlib
from functools import partial
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

#taken from Stack overflow
def md5sum(filename):
    with open(filename, mode='rb') as f:
        d = hashlib.md5()
        for buf in iter(partial(f.read, 128), b''):
            d.update(buf)
    return d.hexdigest()


def test_plink_file_writing():
    """
    This tests byte for byte similarity between plink files that originated from plink.
    and that originated from this package.
    :return:
    """
    print(__file__)
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
    orig_bed_md5sum, orig_bim_md5sum, orig_fam_md5sum = [md5sum(f'{plink_loc}{x}') for x in ['.bed', '.bim', '.fam']]

    #this is where all the new files are made.
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

    # clean up
    os.remove(f"{out_loc}.fam")
    os.remove(f"{out_loc}.bim")
    os.remove(f"{out_loc}.bed")




def test_plink_file_writing_after_variant_prune():
    """
    This tests byte for byte similarity between plink files pruned for variants that originated from plink.
    and that bed fiules pruned for variants originated from this package.
    :return:
    """
    rel_path = '/'.join(('test_resources', 'subset_of_exposure_cohort'))
    if len(__file__.split("/")) > 1:
        plink_loc = "{}/{}".format("/".join(__file__.split("/")[:-1]), rel_path)
    else:
        plink_loc = rel_path

    out_path = '/'.join(('temp_data', 'check_bed_pruning'))
    if len(__file__.split("/")) > 1:
        out_loc = "{}/{}".format("/".join(__file__.split("/")[:-1]), out_path)
    else:
        out_loc = out_path

    plinkfile = utils.PlinkFile(plink_loc)
    plinkfile.read_bed_file_into_numpy_array(allele_2_as_zero=False)


    variants_to_prune = ["rs983287", "rs9789701", "rs9653466", "rs9653467", "rs9677214", "rs962469", "rs9678978",
                         "rs974374", "rs9653468", "rs9973918", "rs972372", "rs9967668", "rs974950", "rs9967745",
                         "rs9679402","rs992153", "rs9677150", "rs956730", "rs997049", "rs9808381", "rs995515",
                         "rs995514", "rs955754", "rs9646944", "rs9679297", "rs953934"]

    variants_to_extract_file = out_loc + '_variants_to_extract'
    with open(variants_to_extract_file, 'w') as f:
        for line in variants_to_prune:
            f.write(f'{line}\n')

    plink_implementation_prepend = f"{out_loc}_plink_implementation"
    subprocess.run(['plink',
                     '--bfile', plink_loc,
                     '--extract', variants_to_extract_file,
                     '--make-bed',
                     '--out', plink_implementation_prepend],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    orig_bed_md5sum, orig_bim_md5sum, orig_fam_md5sum = [md5sum(f'{plink_implementation_prepend}{x}')
                                                         for x in ['.bed', '.bim', '.fam']]


    #do the pruning step in plink
    plinkfile = plinkfile.prune_for_a_list_of_snps(variants_to_prune)

    #this is where all the new files are made.
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

    # clean up
    os.remove(f"{out_loc}.fam")
    os.remove(f"{out_loc}.bim")
    os.remove(f"{out_loc}.bed")

    os.remove(f'{variants_to_extract_file}')

    os.remove(f"{plink_implementation_prepend}.fam")
    os.remove(f"{plink_implementation_prepend}.bim")
    os.remove(f"{plink_implementation_prepend}.bed")
    os.remove(f"{plink_implementation_prepend}.log")
    os.remove(f"{plink_implementation_prepend}.nosex")



def test_plink_file_writing_after_region_prune():
    """
    This tests byte for byte similarity between plink files pruned for variants that originated from plink.
    and that bed files pruned for variants originated from this package.
    :return:
    """
    rel_path = '/'.join(('test_resources', 'subset_of_exposure_cohort'))

    if len(__file__.split("/")) > 1:
        plink_loc = "{}/{}".format("/".join(__file__.split("/")[:-1]), rel_path)
    else:
        plink_loc = rel_path

    out_path = '/'.join(('temp_data', 'check_bed_region_pruning'))
    if len(__file__.split("/")) > 1:
        out_loc = "{}/{}".format("/".join(__file__.split("/")[:-1]), out_path)
    else:
        out_loc = out_path

    plinkfile = utils.PlinkFile(plink_loc)
    plinkfile.read_bed_file_into_numpy_array(allele_2_as_zero=False)

    prune_region = gene_regions.StartEndRegion(['2', 102e6, 103e6])

    plink_implementation_prepend = f"{out_loc}_plink_implementation"
    subprocess.run(['plink',
                     '--bfile', plink_loc,
                     '--chr', '2',
                     '--from-mb', '102',
                     '--to-mb', '103',
                     '--make-bed',
                     '--out', plink_implementation_prepend],
                   check=True,
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)

    orig_bed_md5sum, orig_bim_md5sum, orig_fam_md5sum = [md5sum(f'{plink_implementation_prepend}{x}')
                                                         for x in ['.bed', '.bim', '.fam']]

    #do the pruning step in plink
    plinkfile = plinkfile.prune_for_a_region(prune_region)

    #this is where all the new files are made.
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

    # clean up
    os.remove(f"{out_loc}.fam")
    os.remove(f"{out_loc}.bim")
    os.remove(f"{out_loc}.bed")

    os.remove(f"{plink_implementation_prepend}.fam")
    os.remove(f"{plink_implementation_prepend}.bim")
    os.remove(f"{plink_implementation_prepend}.bed")
    os.remove(f"{plink_implementation_prepend}.log")
    os.remove(f"{plink_implementation_prepend}.nosex")



if __name__ == '__main__':
    test_plink_file_writing()
    test_plink_file_writing_after_variant_prune()
    test_plink_file_writing_after_region_prune()