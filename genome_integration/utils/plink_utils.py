import os
import numpy as np
import sklearn
import subprocess
import bitarray
import copy
from .file_utils import *
from .gcta_ma_utils import *
from .. samples import Sample
from .. variants import SNP
from .. import variants


class BimFile:
    """
    Reads in a bim  file

    Attributes
    ----------

    bim_results: dict
        Values contains SNP information of all variants in the bim file. keys are the SNP name of the bim file.

    bim_results_by_pos: dict
        Values contains SNP information of all variants in the bim file. keys are the pos name "<chr>:<pos>"
        of the bim file.

    snp_names: list
        ordered list of snp names

    Methods
    -------

    add_frq_information(self, file_name):
        Adds frequency information from a plink file (.frq or .freqx)


    """

    def __init__(self, file_name):
        self.bim_results = {}
        self.bim_results_by_pos = {}
        self.snp_names = []
        with open(file_name, 'r') as f:
            for line in f:

                split = [x for x in line[:-1].split() if x != ""]
                tmp = SNP(snp_name=split[1],
                            chromosome=split[0],
                            position=split[3],
                            major_allele=split[5],
                            minor_allele=split[4],
                            minor_allele_frequency=None
                            )

                posname = str(tmp.chromosome) + ":" + str(tmp.position)
                self.snp_names.append(tmp.snp_name)
                self.bim_results[tmp.snp_name] = tmp
                self.bim_results_by_pos[posname] = tmp



    def add_frq_information(self, file_name):
        """
        Adds frequency information from a plink file (.frq or .freqx)

        :param file_name: file name of the plink frq or frqx file.
        :return: self, with added frq information

        """
        with open(file_name, 'r') as f:
            split = f.readline()[:-1].split("\t")

            # frq file.
            if len(split) == 5:
                for line in f:
                    split = [x for x in line.split() if x != ""]
                    snp_name = split[1]
                    try:
                        self.bim_results[snp_name].add_minor_allele_frequency(split[3], split[2], float(split[4]))
                    except KeyError:
                        try:
                            self.bim_results_by_pos[snp_name].add_minor_allele_frequency(split[3], split[2], float(split[4]))
                        except KeyError:
                            continue

            # frqx file
            elif len(split) == 10:
                for line in f:
                    split = [x for x in line.split() if x != ""]
                    snp_name = split[1]
                    a_one_count = int(split[4])*2 + int(split[5])
                    a_two_count = int(split[6])*2 + int(split[5])
                    if a_one_count <= a_two_count:
                        minor = split[2]
                        major = split[3]
                        maf = float(a_one_count) / float(a_one_count + a_two_count)

                    else:
                        minor = split[3]
                        major = split[2]

                        maf = float(a_two_count) / float((a_one_count + a_two_count))
                    try:
                        self.bim_results[snp_name].add_minor_allele_frequency(major, minor, float(maf))
                    except KeyError:
                        try:
                            self.bim_results_by_pos[snp_name].add_minor_allele_frequency(major, minor, float(maf))
                        except:
                            continue
            else:
                RuntimeError("The frq file header was not in any correct formatting.")


class FamSample(Sample):
    """
    Extends the sample class, will also contain family ID and sample ID.

    Attributes
    ----------

    special to this class

    fid: fid from the plink file
    iid: iid from the plink file


    """
    def __init__(self, fid, iid, sex, name, phenotype=None):

        super().__init__(name, phenotype)
        self.fid =  fid
        self.iid = iid
        self.sex = sex



class FamFile:
    """
    implements a fam file

    Attributes
    ----------
    fam_loc: str
        path to the fam file.

    sample_names: list
        sample names of the fam. Currently using <fid>~__~<iid>

    fam_samples: dict
        dictionary with sample names as keys and FamSample as the value.

    """

    def __init__(self, fam_loc):
        self.fam_loc = fam_loc
        self.sample_names = []
        self.fam_samples = {}


        with open(fam_loc, "r") as f:
            for line in f:
                split = line.split()
                sample_name = f'{split[0]}~__~{split[1]}'
                if sample_name in self.sample_names:
                    raise ValueError(f"{self.fam_loc} contains multiple individuals with the same ID.")

                if split[5] == "-9":
                    this_individual = FamSample(split[0], split[1], split[4], sample_name)
                else:
                    this_individual = FamSample(split[0], split[1], split[4], sample_name, split[5])

                self.sample_names.append(this_individual.name)
                self.fam_samples[this_individual.name] = this_individual




class PlinkFile:
    """
    A reader for the plink bed file format.

    requires a valid location of the bed file.

    Attributes
    ----------
    bed_loc: str
        of the bed file location

    bim_loc: str
        of the bim file location

    fam_loc: str
        of the fam file location

    genotypes: float numpy array of variants on columns samples on the rows.
        Numpy array of genotypes default is number of minor alleles per variant and position.
        not initialized, need to call read_bed_file_into_numpy_array

    bim_data: BimFile object.
        contains all the variants from the bim file.

    fam_data: FamFile object
        contains all the Samples from the bim file

    _decoder: dict
        bit dict of the file encoding.

    Methods
    -------
    read_bed_file_into_numpy_array(self, allele_2_as_zero=True, missing_encoding=3)
        reads in the bed file takes about 5 seconds for 5,000 individuals ~14,000 variants.
        saves it to the class.

    """

    def __init__(self, bfile_location):

        self.bed_loc = bfile_location + ".bed"
        self.bim_loc = bfile_location + ".bim"
        self.fam_loc = bfile_location + ".fam"

        self.genotypes = None

        if not os.path.exists(self.bed_loc):
            raise FileNotFoundError(f"Could not find bimfile: {self.bed_loc}")

        if not os.path.exists(self.bim_loc):
            raise FileNotFoundError(f"Could not find bimfile: {self.bim_loc}")

        if not os.path.exists(self.bim_loc):
            raise FileNotFoundError(f"Could not find bimfile: {self.fam_loc}")

        #variants
        self.bim_data = BimFile(self.bim_loc)

        if len(self.bim_data.snp_names) > 100000:
            print("Warning, Plinkfile with more than 100,000 variants will consume a lot of time and memory.")

        self.fam_data = FamFile(self.fam_loc)
        self._decoder = {0: bitarray.bitarray('00'),  #homoz A1 (usually minor)
                         1: bitarray.bitarray('01'),  #heteroz
                         2: bitarray.bitarray('11'),  #homoz A2 (usually major)
                         3: bitarray.bitarray('10'),  # missing
                         }

    def read_bed_file_into_numpy_array(self, allele_2_as_zero=True, missing_encoding=3):
        """
        Reads a bed file into a numpy array

        :param allele_2_as_zero: bool
            allele 2 (often the major allele in plink) is encoded as zero, making an increase in the minor allele an
            increase in number. This is opposite to the plinkio encoding.
            But this makes the minor allele often also the effect allele.

        :param missing_encoding: int
            encodes missing values as three, plinkio default.
        :return: n indivduals by m variants numpy array (floats) of genotypes.
        """
        bed_file = open(self.bed_loc, "rb")
        magick = bed_file.read(3)
        if magick != b'l\x1b\x01':
            raise ValueError("Plink file magic string is not correct.")

        all_bytes = bed_file.read()

        num_variants = len(self.bim_data.snp_names)
        num_individuals = len(self.fam_data.sample_names)

        bytes_to_read = int(np.ceil(num_individuals / 4))

        if len(all_bytes) != (bytes_to_read*num_variants):
            ValueError(f"{self.bed_loc} has an incorrect number of bytes: expected {(bytes_to_read*num_variants)}, found: {len(all_bytes)}")
        bed_file.close()

        genotypes = np.zeros( (num_individuals, num_variants), dtype=np.uint8)

        offset = 0
        for i in range(num_variants):

            bits = bitarray.bitarray(endian="little")
            bits.frombytes(all_bytes[offset:(offset+bytes_to_read)])
            array = bits.decode(self._decoder)[:num_individuals]
            genotypes[:,i] = array
            offset+=bytes_to_read

        #convert genotypes if necessary (default)
        if allele_2_as_zero:
            tmp_geno = copy.deepcopy(genotypes)
            tmp_geno[genotypes == 2] = 0
            tmp_geno[genotypes == 0] = 2
            genotypes = tmp_geno

        genotypes = np.array(genotypes, dtype=float)
        genotypes[genotypes == 3] = missing_encoding

        self.genotypes = genotypes
        return genotypes





def plink_isolate_clump(bed_file, associations, threshold, r_sq=0.5  ,tmploc="", return_snp_file=False):
    """
    will prune for ld in a list of snps. from a bed file location.
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

    subprocess.run(["plink --bfile " + bed_file +
                    " --clump " + clump_file +
                    " --clump-p1 " + str(threshold) +
                    " --clump-r2 " + str(r_sq) +
                    " --clump-kb 1000 " +
                    " --out " + tmp_plink_out],
                   shell=True,
                   check=True,
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)

    clumpdat = read_newline_separated_file_into_list(tmp_plink_out + ".clumped")

    snps_to_keep = [x.split()[2] for x in clumpdat[1:] if len(x.split())]

    if return_snp_file:
        write_list_to_newline_separated_file(snps_to_keep, snp_out)

    subprocess.run(["rm " + tmp_plink_out + ".* " + clump_file ], shell=True, check=True)

    if return_snp_file:
        return snps_to_keep, snp_out

    else:
        return snps_to_keep


def isolate_snps_of_interest_make_bed(ma_file, exposure_name, b_file,
                                      snp_file_out, plink_files_out, calculate_ld = False):
    """

    Isolate snps of interest for a gene, and make a bed file

    :param ma_file:
    :param exposure_name:
    :param b_file:
    :param snp_file_out:
    :param plink_files_out:
    :param calculate_ld:
    :return: the name_of the bedfile with only the snps
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
                             stdout=subprocess.DEVNULL,  # to DEVNULL, because plink saves a log of everything
                             stderr=subprocess.DEVNULL
                             )

        bim_file = BimFile(plink_files_out + '.bim')

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
                             stdout=subprocess.DEVNULL, # to DEVNULL, because plink saves a log of everything
                             stderr = subprocess.DEVNULL
                             )

        bim_file = BimFile(plink_files_out + '.bim')

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



def score_and_assess_auc(genetic_associations, bed_file, tmp_file = "tmp_score", p_value_thresh = 1.0, resolution = 500):
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
