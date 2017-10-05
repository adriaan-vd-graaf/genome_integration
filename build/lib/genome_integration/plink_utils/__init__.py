"""
These functions are used to use plink that should be installed in your system.
"""

import subprocess
from .. import gcta_utils
from .. import file_utils
from .. import variants


def plink_isolate_clump(bed_file, associations, threshold, r_sq=0.5  ,tmploc="", return_snp_file=False):
    """

    will prune for ld in a list of snps.
    will output a list after prune.

    :param bed_file:
    :param associations:
    :return: list of snps after prune
    """
    clump_file = tmploc +  "_" + str(threshold) +"_clump_file.txt"
    tmp_plink_out = tmploc +  "_" + str(threshold) + "_plink_out"
    snp_out = tmploc + "_" + str(threshold) + "clumped.txt"

    association_lines = ["SNP\tP"]
    [association_lines.append("{}\t{}".format(
        associations[x].snp_name,
        associations[x].wald_p_val
    )) for x in associations.keys()]

    file_utils.write_list_to_newline_separated_file(association_lines, clump_file)

    subprocess.run(["plink --bfile " + bed_file  +
                    " --clump " + clump_file +
                    " --clump-p1 " + str(threshold) +
                    " --clump-r2 " + str(r_sq) +
                    " --clump-kb 1000 "+
                    " --out " + tmp_plink_out], shell=True, check=True, stdout=subprocess.DEVNULL)

    clumpdat = file_utils.read_newline_separated_file_into_list(tmp_plink_out + ".clumped")

    snps_to_keep = [x.split()[2] for x in clumpdat[1:] if len(x.split())]

    if return_snp_file:
        file_utils.write_list_to_newline_separated_file(snps_to_keep, snp_out)

    subprocess.run(["rm " + tmp_plink_out + ".* " + clump_file ], shell=True, check=True)

    if return_snp_file:
        return snps_to_keep, snp_out

    else:
        return snps_to_keep


def make_ld_mat_from_genetic_associations(genetic_associations, bfile, tmp_dir, thread_num=1):

    bfile_out = tmp_dir + "_bed_file"
    snp_out = tmp_dir + "snps"

    file_utils.write_list_to_newline_separated_file([genetic_associations[x].snp_name for x in genetic_associations.keys()], snp_out)

    tmp = subprocess.run(['plink',
                          '--bfile', bfile,
                          '--extract', snp_out,
                          '--make-bed',
                          '--r', 'square',
                          '--threads', str(int(thread_num)),
                          '--out', bfile_out
                          ],
                         check=True,
                         stdout=subprocess.DEVNULL  # to DEVNULL, because plink saves a log of everything
                         )

    bim_file = variants.BimFile(bfile_out+ '.bim')
    bim_file.add_ld_mat(bfile_out + '.ld')

    subprocess.run('rm -f {} {}*'.format(snp_out, bfile_out), shell=True)

    return bim_file



def isolate_snps_of_interest_make_bed(ma_file, exposure_name, b_file, snp_file_out, plink_files_out, calculate_ld = False):
    """
    Isolate snps of interest for a gene, and make a bed file

    :return:
    the name_of the bedfile with only the snps
    """

    ma_data = gcta_utils.MaFile(ma_file, exposure_name)

    # write the snps to isolate
    file_utils.write_list_to_newline_separated_file(ma_data.snp_names(no_palindromic=True), snp_file_out)
    if calculate_ld:
        # now run plink to isolate the files, and return the snplist, plink filename and eqtl ma file.
        tmp = subprocess.run(['plink',
                              '--bfile', b_file,
                              '--extract', snp_file_out,
                              '--make-bed',
                              '--r', 'square',
                              '--out', plink_files_out
                             ],
                             check=True,
                             stdout=subprocess.DEVNULL  # to DEVNULL, because plink saves a log of everything
                             )

        bim_file = variants.BimFile(plink_files_out + '.bim')
        bim_file.add_ld_mat(plink_files_out + '.ld')

        return ma_data, bim_file

    else:
        # now run plink to isolate the files, and return the snplist, plink filename and eqtl ma file.
        tmp = subprocess.run(['plink',
                              '--bfile', b_file,
                              '--extract', snp_file_out,
                              '--make-bed',
                              '--out', plink_files_out
                              ],
                             check=True,
                             stdout=subprocess.DEVNULL  # to DEVNULL, because plink saves a log of everything
                             )

        bim_file = variants.BimFile(plink_files_out + '.bim')

        return ma_data, bim_file