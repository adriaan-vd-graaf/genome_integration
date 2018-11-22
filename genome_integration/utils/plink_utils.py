
import numpy as np
import sklearn
import subprocess
from . import *
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

    write_list_to_newline_separated_file(association_lines, clump_file)

    subprocess.run(["plink --bfile " + bed_file  +
                    " --clump " + clump_file +
                    " --clump-p1 " + str(threshold) +
                    " --clump-r2 " + str(r_sq) +
                    " --clump-kb 1000 "+
                    " --out " + tmp_plink_out], shell=True, check=True, stdout=subprocess.DEVNULL)

    clumpdat = read_newline_separated_file_into_list(tmp_plink_out + ".clumped")

    snps_to_keep = [x.split()[2] for x in clumpdat[1:] if len(x.split())]

    if return_snp_file:
        write_list_to_newline_separated_file(snps_to_keep, snp_out)

    subprocess.run(["rm " + tmp_plink_out + ".* " + clump_file ], shell=True, check=True)

    if return_snp_file:
        return snps_to_keep, snp_out

    else:
        return snps_to_keep


def make_ld_mat_from_genetic_associations(genetic_associations, bfile, tmp_dir, thread_num=1):

    bfile_out = tmp_dir + "_bed_file"
    snp_out = tmp_dir + "snps"

    write_list_to_newline_separated_file([genetic_associations[x].snp_name for x in genetic_associations.keys()], snp_out)

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

    ma_data = MaFile(ma_file, exposure_name)

    # write the snps to isolate
    write_list_to_newline_separated_file(ma_data.snp_names(no_palindromic=True), snp_file_out)
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


def score_individuals(genetic_associations, bed_file, tmp_file = "tmp_score", p_value_thresh = 1):
    """
    Used to score individual.
    :param genetic_associations:
    :param bed_file: prepend of a bed file
    :param tmp_file: prepend of temporary files.
    :param p_value_thresh: p value threshold of which the genetic associations should be part of.
    :return: dict with keys corresponding to individuals,
            values: tuple with the phenotype [0]  and score [1]  of the individual.
    """


    file_for_scoring = tmp_file + "_snps_beta.txt"
    pos_name_scoring = tmp_file + "_posname_beta.txt"
    prepend_for_plink = tmp_file + "_score"

    with open(file_for_scoring, "w") as f:
        for snp in genetic_associations.keys():
            tmp_assoc = genetic_associations[snp]
            if tmp_assoc.wald_p_val < p_value_thresh:
                f.write("{}\t{}\t{}\n".format(tmp_assoc.snp_name, tmp_assoc.minor_allele, tmp_assoc.beta))

    with open(pos_name_scoring, "w") as f:
        for snp in genetic_associations.keys():
            tmp_assoc = genetic_associations[snp]
            if tmp_assoc.wald_p_val < p_value_thresh:
                f.write("{}\t{}\t{}\n".format("{}:{}".format(tmp_assoc.chromosome, tmp_assoc.position), tmp_assoc.minor_allele, tmp_assoc.beta))
    try:
        subprocess.run(["plink",
                        "--allow-no-sex",
                        "--bfile", bed_file,
                        "--score",  file_for_scoring,
                        "--out", prepend_for_plink + ".snp_name"
                        ],
                       check=True,
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL
                       )
        profile_loc = prepend_for_plink + ".snp_name.profile"

    except subprocess.CalledProcessError:
        # something went wrong. Now trying it with snps which have their name as position.
        subprocess.run(["plink",
                        "--allow-no-sex",
                        "--bfile", bed_file,
                        "--score", pos_name_scoring,
                        "--out", prepend_for_plink + ".pos_name"
                        ],
                       check=True,
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL
                       )
        profile_loc = prepend_for_plink + ".pos_name.profile"

    # scoring done, now read the file.
    pheno_score = {}
    with open(profile_loc, "r") as f:
        f.readline()
        for line in f:
            split = line.split()
            pheno_score[split[1]] = (float(split[2]), float(split[5]))

    subprocess.run(['rm -f ' + file_for_scoring + " " + pos_name_scoring + " " + prepend_for_plink + ".*"], shell=True, check=True)

    return pheno_score



def score_and_assess_auc(genetic_associations, bed_file, tmp_file = "tmp_score", p_value_thresh = 1, resolution = 500):
    """
    Using scoring, we determine the auc

    :param genetic_associations:
    :param bed_file:
    :param tmp_file:
    :param p_value_thresh:
    :param resolution:
    :return:
    """

    pheno_score = score_individuals(genetic_associations, bed_file, tmp_file, p_value_thresh)

    pheno = np.array([pheno_score[x][0] for x in pheno_score.keys()])
    scores = np.array([pheno_score[x][1] for x in pheno_score.keys()])

    thresholds = np.arange(min(scores), max(scores), (max(scores) - min(scores)) / resolution)

    # Using 2 as the phenotype celiac disease.
    tpr = [sum(pheno[scores > x] == 2.0) / sum(pheno == 2.0) for x in thresholds]
    fpr  = [sum(pheno[scores > x] == 1.0) / sum(pheno == 1.0) for x in thresholds]

    auc = sklearn.metrics.auc(fpr,tpr)


    return tpr, fpr, auc
