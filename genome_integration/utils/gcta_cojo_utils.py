import subprocess
from . import *

from .plink_utils import *



class CojoCmaFile:
    """
    Contains the GCTA COJO CMA file results.
    The conditional joint associations are located in the ma results part of the file.

    Attributes
    ----------

    name: str
        name of the association

    ma_results: dict
        dict with snp_names as keys, and the CojoCmaLine as values.

    """

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


class CojoCmaLine(association.GeneticAssociation):
    """
    This is a helper class arounc the GeneticAssociation class.
    but adds the information that GCTA also keeps.

    From the PCTG, GCTA documentation:
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

    Attributes specific to this class, rest is inherited from GeneticAssociation:

    beta_initial: float
        initial beta from input
    se_intitial: float
        initial se from input
    p_intitial: float
        initial p value from input
    freq_geno: float
        frequency of the reference populaiton
    n_estimated:
        number of estimated individuals in the population


    """


    def __init__(self, line, name):
        split = [x for x in line.split() if x != ""]
        if float(split[4]) > 0.5:
            major = split[3]
            minor = None  #if I don't know it, it could be N, maybe change later. Use the add_snp_data to update with info
            beta = -1 * float(split[10])
            frq = 1 - float(split[3])
        else:
            minor = split[3]
            major = None
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



def do_gcta_cojo_slct(bfile_prepend, ma_file, out_prepend, p_val='1e-8', maf='0.01', gc=1.0):
    """
    Doeas GCTA COJO stepwise selection (no joint effects)


    :param bfile_prepend: bedfile location
    :param ma_file: ma file location
    :param out_prepend: where to output
    :param p_val: p value threshold for stepwise selection
    :param maf: minor allele frequency threshold for stepwise selection
    :param gc: genomic correction factor, default is 1.0. Make sure to check this in your associations
    :return: CojoCmaFile object with the results.
    """
    base_file = None
    for chr in range(1, 23):
        std_out = open(out_prepend + '.out', 'w')
        std_err = open(out_prepend + '.err', 'w')
        subprocess.run(['gcta64',
                        '--chr', str(chr),
                        '--bfile', bfile_prepend,
                        '--cojo-file', ma_file,
                        '--cojo-slct',
                        '--cojo-gc', str(gc),
                        '--out', out_prepend,
                        '--cojo-p', p_val,
                        '--maf', maf,
                        '--thread-num', '1'
                        ],
                       stdout=std_out,
                       stderr=std_err,
                       check=True
                       )
        std_out.close()
        std_err.close()
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






def do_gcta_cojo_joint(bfile_prepend, ma_file, out_prepend, p_val='1e-8', maf='0.01', gc=1.0):
    """
    Does GCTA COJO stepwise selection and joint effects


    :param bfile_prepend: bedfile location
    :param ma_file: ma file location
    :param out_prepend: where to output
    :param p_val: p value threshold for stepwise selection
    :param maf: minor allele frequency threshold for stepwise selection
    :param gc: genomic correction factor, default is 1.0. Make sure to check this in your associations
    :return: CojoCmaFile object with the results.
    """

    std_out = open(out_prepend + '.out', 'w')
    subprocess.run(['gcta64',
                    '--bfile', bfile_prepend,
                    '--cojo-file', ma_file,
                    '--cojo-slct',
                    '--out', out_prepend,
                    '--cojo-gc', str(gc),
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
                                         p_val_thresh=0.05, maf=0.01, calculate_ld = False,
                                         clump=False, create_tmp_subset_of_bed=True):
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
        clumped_snps = plink_isolate_clump(bed_file=bfile,
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


    ma_lines = [make_gcta_ma_header()]

    [
        ma_lines.append(make_gcta_ma_line(genetic_associations[x])) for x in genetic_associations
        if genetic_associations[x].snp_name in clumped_snps
    ]

    write_list_to_newline_separated_file(ma_lines, ma_name)
    if create_tmp_subset_of_bed:
        try:
            isolate_snps_of_interest_make_bed(ma_file=ma_name, exposure_name=gene_name, b_file=bfile,
                                                          snp_file_out=snp_out, plink_files_out=plink_pruned,
                                                          calculate_ld=calculate_ld)

        except Exception as x:
            print("isolating snps raised an exception while processing " + gene_name )
            # subprocess.run(["rm -f {} {} {}*".format(ma_name, snp_out, plink_pruned)], shell=True, check=True)
            raise x
    else:
        plink_pruned = bfile


    try:
        cojo_eqtl = do_gcta_cojo_slct(plink_pruned, ma_name, cojo_out, p_val='{:6.2e}'.format(p_val_thresh), maf='{:8.6f}'.format(maf))
    except Exception as x:
        print("GCTA cojo raised an exceptqion while processing" + gene_name)
        # subprocess.run(["rm {} {} {}* {}*".format(ma_name, snp_out, plink_pruned, cojo_out)], shell=True, check=True)
        raise x
    if create_tmp_subset_of_bed:
        subprocess.run(["rm -f {} {} {}* {}*".format(ma_name,snp_out,plink_pruned,cojo_out)], shell=True, check = True)
    else:
        subprocess.run(["rm -f {} {} {}*".format(ma_name, snp_out, cojo_out)], shell=True, check=True)
    return cojo_eqtl


def do_gcta_cojo_joint_on_genetic_associations(genetic_associations, bfile, tmp_prepend,
                                         p_val_thresh=0.05,
                                         maf=0.01,
                                         calculate_ld = False,
                                         clump=False,
                                         _keep_ma_files=False):
    """
    :param genetic_associations: a dict of genetic associations,  keys should be explantory name
    :param bfile: plink bed file
    :param tmp_prepend: temporary name of files where to store.
    :param p_val_thresh: p value threshold as a float
    :param maf: minor allele frequency as a float
    :param _keep_ma_files: This is used for testing. MA files are used to know the exact floating point input.
    :return: Cojo results a Cojo CMA file object, which is an extension of the geneticassociation file.
    """


    # define the names.
    ma_name = tmp_prepend + "_temp.ma"
    snp_out = tmp_prepend + "_snp_list.txt"
    plink_pruned = tmp_prepend + "_plink_pruned"
    cojo_out = tmp_prepend + "_cojo_out"

    #clump.
    if clump:
        clumped_snps = plink_isolate_clump(bed_file=bfile,
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


    ma_lines = [make_gcta_ma_header()]

    [
        ma_lines.append(make_gcta_ma_line(genetic_associations[x]))
        for x in genetic_associations
        if genetic_associations[x].snp_name in clumped_snps
    ]

    write_list_to_newline_separated_file(ma_lines, ma_name)

    try:
        isolate_snps_of_interest_make_bed(ma_file=ma_name, exposure_name=gene_name, b_file=bfile,
                                                      snp_file_out=snp_out, plink_files_out=plink_pruned,
                                                      calculate_ld=calculate_ld)

    except Exception as x:
        print("isolating snps raised an exception while processing " + gene_name )
        subprocess.run(["rm -f {} {} {}*".format(ma_name, snp_out, plink_pruned)], shell=True, check=True)
        raise x


    try:
        cojo_eqtl = do_gcta_cojo_joint(plink_pruned, ma_name, cojo_out, p_val='{:6.2e}'.format(p_val_thresh), maf='{:8.6f}'.format(maf))
    except Exception as x:
        print("GCTA cojo raised an exception while processing " + gene_name)
        subprocess.run(["rm {} {} {}* {}*".format(ma_name, snp_out, plink_pruned, cojo_out)], shell=True, check=True)
        raise x

    if _keep_ma_files:
        subprocess.run(["rm -f {} {}* {}*".format(snp_out,plink_pruned,cojo_out)], shell=True, check = True)
    else:
        subprocess.run(["rm -f {} {} {}* {}*".format(ma_name,snp_out,plink_pruned,cojo_out)], shell=True, check = True)

    return cojo_eqtl


def make_gcta_ma_header():
    """
    Will create an ma header. for GCTA-COJO

    :return: String with an ma file header.
    """
    return "SNP\tA1\tA2\tfreq\tb\tse\tp\tN"


def make_gcta_ma_line(genetic_association):
    """

    Makes a GCTA line of the genetic variant.

    Will only return a string not newline ended,
    will not write to a file, the user is expected to do this himself.

    :param genetic association class object.
    :return tab separated string that can be part of ma file:
    """

    # make sure the data that we need is available.
    if not genetic_association.has_position_data or not genetic_association.has_allele_data or not genetic_association.has_frequency_data:
        raise RuntimeError("Cannot write an Ma line. Does not contain the necessary data")

    if genetic_association.wald_p_val == None:
        raise RuntimeError("No p value present")

    return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(genetic_association.snp_name,
                                                   genetic_association.minor_allele,
                                                   genetic_association.major_allele,
                                                   genetic_association.minor_allele_frequency,
                                                   genetic_association.beta,
                                                   genetic_association.se,
                                                   genetic_association.wald_p_val,
                                                   genetic_association.n_observations)