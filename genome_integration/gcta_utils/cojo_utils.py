import re
import time
import subprocess
import numpy as np
from .. import variants
from .. import association
from .. import file_utils
from .. import plink_utils
from collections import Counter

class CojoCmaFile:
    def __init__(self, file_loc, name):
        self.name = name
        self.ma_results = {}
        with open(file_loc, 'r') as f:
            f.readline()
            for line in f:
                try:
                    tmp = CojoCmaLine(line, name)
                    self.ma_results[tmp.snp_name] = tmp
                except ValueError as x:
                    # print("Could not find a valid line in the following line: {}".format(line))
                    continue

    def snps_with_data(self):
        return set(self.ma_results.keys())

    def write_cojo_result(self, file_name):
        lines = [''] * (len(self.ma_results.keys())+1)
        lines[0] = "snp_name\tchr\tbp\tbeta\tse\tp_val\tassoc_name"
        indice = 1
        for i in list(self.ma_results.keys()):
            tmp = self.ma_results[i]
            lines[indice] = '\t'.join([tmp.snp_name,
                               str(tmp.chromosome),
                               str(tmp.position),
                               str(tmp.beta),
                               str(tmp.se),
                               str(tmp.wald_p_val),
                               str(self.name)
                               ]
                              )
            indice += 1


        file_utils.write_list_to_newline_separated_file(lines, file_name)


class CojoCmaLine(association.GeneticAssociation):
    """
    From the GCTA website:
    Columns are:

    chromosome;
    SNP;
    physical position;
    frequency of the effect allele in the original data;
    the effect allele;
    effect size,
    standard error and
    p-value from the original GWAS or meta-analysis;
    estimated effective sample size;
    frequency of the effect allele in the reference sample;
    effect size,
    standard error and
    p-value from a joint analysis of all the selected SNPs;
    LD correlation between the SNP i and SNP i + 1 for the SNPs on the list.

    """

    __slots__ = ['snp_name', 'chromosome', 'position', 'major_allele', 'minor_allele', 'minor_allele_frequency',
                 'has_position_data', 'has_allele_data', 'has_frequency_data', 'dependent_name', 'explanatory_name',
                 'beta', 'se', 'n_observations', 'r_squared', 'z_score', 'wald_p_val', 'snp', 'beta_initial',
                 'se_initial', 'p_initial', 'freq_geno', 'n_estimated']

    def __init__(self, line, name):
        split = [x for x in line.split() if x != ""]
        if float(split[4]) > 0.5:
            major = split[3]
            minor = "N" #if I don't know it, it could be N, maybe change later. Use the add_snp_data to update with info
            beta = -1 * float(split[10])
            frq = 1 - float(split[3])
        else:
            minor = split[3]
            major = "N"
            beta = float(split[10])
            frq = float(split[4])


        super().__init__(
            dependent_name=name,
            explanatory_name=split[1],
            n_observations=float(split[9]),
            beta=beta,
            se=float(split[11]),
            r_squared=None,
            chromosome=split[0],
            position=int(split[2]),
            major_allele=major,
            minor_allele=minor,
            minor_allele_frequency=frq
        )

        self.beta_initial = float(split[5])
        self.se_initial = float(split[6])
        self.p_initial = float(split[7]) # initial p value.
        self.freq_geno  = float(split[8]) #frequency from reference population.
        self.n_estimated = float(split[9])

        self.wald_p_val = float(split[12])


class CojoLdrFile:

    def __init__(self, file_loc, name):
        self.name = name
        with open(file_loc, 'r') as f:
            # first line contains the SNPs
            self.snps = [variants.BaseSNP(x) for x in f.readline()[:-1].split() if (x != '') and (x != 'SNP')]
            self.snp_names = [x.snp_name for x in self.snps]
            self.ld_mat = np.zeros((len(self.snps),len(self.snps)))
            indice = 0
            for line in f:
                self.ld_mat[:,indice] = [float(x) for x in line.split()[1:] if x != '']
                indice += 1

    def get_ld_mat(self):
        return self.ld_mat, self.snp_names

    # this may need to be rewritten, to ensure the snps go from lowest bp position to highest.
    # Not checked right now

    #right now assuming there are no trans effects


    def write_ld_mat_gg(self, filename, bim_data):
        snpnames = [x.snp_name for x in self.snps]

        position = []

        for i in snpnames:
            try:
                position.append(bim_data.bim_results[i].position)
            except:
                position.append(np.NaN)

        ordering = np.argsort(np.array(position))
        string_list = ['chr\tbp\tsnp_name\t' + '\t'.join(np.array(snpnames)[ordering])]
        for i in ordering:
            try:
                tmp = bim_data.bim_results[snpnames[i]].chromosome + '\t' \
                      + bim_data.bim_results[snpnames[i]].position + '\t' \
                      + snpnames[i] + '\t'
            except:
                tmp = 'NA\tNA\t' + snpnames[i] + '\t'

            tmp += '\t'.join([str(x) for x in self.ld_mat[i, ordering]])
            string_list.append(tmp)

        file_utils.write_list_to_newline_separated_file(string_list, filename)


# todo make this work into a single function that accepts a dict and a temporary folder.
def do_gcta_cojo_slct(bfile_prepend, ma_file, out_prepend, p_val='1e-8', maf='0.01'):

    base_file = None
    for chr in range(1, 23):
        std_out = open(out_prepend + '.out', 'w')
        subprocess.run(['gcta64',
                        '--chr', str(chr),
                        '--bfile', bfile_prepend,
                        '--cojo-file', ma_file,
                        '--cojo-slct',
                        '--out', out_prepend,
                        '--cojo-p', p_val,
                        '--maf', maf,
                        '--thread-num', '1'
                        ],
                       stdout=std_out,
                       check=True
                       )
        std_out.close()
        try:
            tmp_cojo = CojoCmaFile(out_prepend + ".jma.cojo", out_prepend)
        except:
            continue

        if base_file is None:
            base_file = tmp_cojo
        else:
            for i in tmp_cojo.ma_results.keys():
                base_file.ma_results[i] = tmp_cojo.ma_results[i]

    return base_file






def do_gcta_cojo_joint(bfile_prepend, ma_file, out_prepend, p_val='1e-8', maf='0.01'):


    std_out = open(out_prepend + '.out', 'w')
    subprocess.run(['gcta64',
                    '--bfile', bfile_prepend,
                    '--cojo-file', ma_file,
                    '--cojo-slct',
                    '--out', out_prepend,
                    '--cojo-p', p_val,
                    '--maf', maf,
                    '--thread-num', '1'
                    ],
                   stdout=std_out,
                   check=True
                   )
    std_out.close()

    tmp_cojo = CojoCmaFile(out_prepend + ".jma.cojo", out_prepend)

    return tmp_cojo


def do_gcta_cojo_on_genetic_associations(genetic_associations, bfile, tmp_prepend,
                                         p_val_thresh=0.05, maf=0.01, calculate_ld = False, clump=False):
    """
    :param genetic_associations: a dict of genetic associations,  keys should be explantory name
    :param bfile: plink bed file
    :param tmp_prepend: temporary name of files where to store.
    :param p_val_thresh: p value threshold as a float
    :param maf: minor allele frequency as a float

    :return: Cojo results a Cojo CMA file object, which is an extension of the geneticassociation file.
    """


    # define the names.
    ma_name = tmp_prepend + "_temp.ma"
    snp_out = tmp_prepend + "_snp_list.txt"
    plink_pruned = tmp_prepend + "_plink_pruned"
    cojo_out = tmp_prepend + "_cojo_out"

    #clump.
    if clump:
        clumped_snps = plink_utils.plink_isolate_clump(bed_file=bfile,
                                                   associations=genetic_associations,
                                                   threshold=p_val_thresh,
                                                   r_sq=0.5,
                                                   tmploc=tmp_prepend
                                                   )
    else:
        clumped_snps = list(genetic_associations.keys())

    snps = list(genetic_associations.keys())
    try:
        gene_name = genetic_associations[snps[0]].dependent_name.decode("utf8")
    except AttributeError:
        gene_name = genetic_associations[snps[0]].dependent_name


    ma_lines = [genetic_associations[snps[0]].make_gcta_ma_header()]

    [
        ma_lines.append(genetic_associations[x].make_gcta_ma_line()) for x in genetic_associations
        if genetic_associations[x].explanatory_name in clumped_snps
    ]

    file_utils.write_list_to_newline_separated_file(ma_lines, ma_name)

    try:
        plink_utils.isolate_snps_of_interest_make_bed(ma_file=ma_name, exposure_name=gene_name, b_file=bfile,
                                                      snp_file_out=snp_out, plink_files_out=plink_pruned,
                                                      calculate_ld=calculate_ld)

    except Exception as x:
        print("isolating snps raised an exception while processing " + gene_name )
        # subprocess.run(["rm -f {} {} {}*".format(ma_name, snp_out, plink_pruned)], shell=True, check=True)
        raise x


    try:
        cojo_eqtl = do_gcta_cojo_slct(plink_pruned, ma_name, cojo_out, p_val='{:6.2e}'.format(p_val_thresh), maf='{:8.6f}'.format(maf))
    except Exception as x:
        print("GCTA cojo raised an exceptqion while processing" + gene_name)
        # subprocess.run(["rm {} {} {}* {}*".format(ma_name, snp_out, plink_pruned, cojo_out)], shell=True, check=True)
        raise x

    subprocess.run(["rm -f {} {} {}* {}*".format(ma_name,snp_out,plink_pruned,cojo_out)], shell=True, check = True)

    return cojo_eqtl


def do_gcta_cojo_joint_on_genetic_associations(genetic_associations, bfile, tmp_prepend,
                                         p_val_thresh=0.05, maf=0.01, calculate_ld = False, clump=False):
    """
    :param genetic_associations: a dict of genetic associations,  keys should be explantory name
    :param bfile: plink bed file
    :param tmp_prepend: temporary name of files where to store.
    :param p_val_thresh: p value threshold as a float
    :param maf: minor allele frequency as a float

    :return: Cojo results a Cojo CMA file object, which is an extension of the geneticassociation file.
    """


    # define the names.
    ma_name = tmp_prepend + "_temp.ma"
    snp_out = tmp_prepend + "_snp_list.txt"
    plink_pruned = tmp_prepend + "_plink_pruned"
    cojo_out = tmp_prepend + "_cojo_out"

    #clump.
    if clump:
        clumped_snps = plink_utils.plink_isolate_clump(bed_file=bfile,
                                                   associations=genetic_associations,
                                                   threshold=p_val_thresh,
                                                   r_sq=0.5,
                                                   tmploc=tmp_prepend
                                                   )
    else:
        clumped_snps = list(genetic_associations.keys())


    snps = list(genetic_associations.keys())
    try:
        gene_name = genetic_associations[snps[0]].dependent_name.decode("utf8")
    except AttributeError:
        gene_name = genetic_associations[snps[0]].dependent_name


    ma_lines = [genetic_associations[snps[0]].make_gcta_ma_header()]

    [
        ma_lines.append(genetic_associations[x].make_gcta_ma_line())
        for x in genetic_associations
        if genetic_associations[x].explanatory_name in clumped_snps
    ]

    file_utils.write_list_to_newline_separated_file(ma_lines, ma_name)

    try:
        plink_utils.isolate_snps_of_interest_make_bed(ma_file=ma_name, exposure_name=gene_name, b_file=bfile,
                                                      snp_file_out=snp_out, plink_files_out=plink_pruned,
                                                      calculate_ld=calculate_ld)

    except Exception as x:
        print("isolating snps raised an exception while processing " + gene_name )
        # subprocess.run(["rm -f {} {} {}*".format(ma_name, snp_out, plink_pruned)], shell=True, check=True)
        raise x


    try:
        cojo_eqtl = do_gcta_cojo_joint(plink_pruned, ma_name, cojo_out, p_val='{:6.2e}'.format(p_val_thresh), maf='{:8.6f}'.format(maf))
    except Exception as x:
        print("GCTA cojo raised an exceptqion while processing " + gene_name)
        # subprocess.run(["rm {} {} {}* {}*".format(ma_name, snp_out, plink_pruned, cojo_out)], shell=True, check=True)
        raise x

    subprocess.run(["rm -f {} {} {}* {}*".format(ma_name,snp_out,plink_pruned,cojo_out)], shell=True, check = True)

    return cojo_eqtl

