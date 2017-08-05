import re
import time
import subprocess
import numpy as np
from .. import variants
from .. import association
from .. import file_utils


class CojoCmaFile:
    def __init__(self, file_loc, name):
        self.name = name
        self.ma_results = {}
        with open(file_loc, 'r') as f:
            f.readline()
            for line in f:
                tmp = CojoCmaLine(line, name)
                self.ma_results[tmp.snp_name] = tmp

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

    ##make sure the process is finished.
    regex = re.compile("^Analysis finished:.*")

    #make sure the log file is valid, assuming the other files are valid as well.
    for i in range(2):
        log_file = file_utils.read_newline_separated_file_into_list(out_prepend + '.out')
        if sum([regex.match(x) != None for x in log_file]) != 1:
            time.sleep(1)
        else:
            break
        if i == 9:
            raise IOError('Cojo analysis was finished, but did not find valid log file here:' + out_prepend + '.out')

    return CojoCmaFile(out_prepend + ".jma.cojo", out_prepend)
