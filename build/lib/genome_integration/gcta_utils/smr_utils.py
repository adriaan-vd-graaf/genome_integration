import scipy.stats

class SmrMaFile:
    def __init__(self, file_loc, name):

        self.name = name
        self.ma_results = {}
        with open(file_loc, 'r') as f:
            f.readline()
            for line in f:
                tmp = SmrResult(line)
                self.ma_results[tmp.probe_id] = tmp

    def snp_names(self):
        return(self.ma_results.keys())


class SmrResult:
    def __init__(self, line=None, outcome_result=None, exposure_result=None):
        if line != None:
            split = [x for x in line[:-1].split() if x != ""]

            self.probe_id = split[0]
            self.probe_chr = int(split[1])
            self.gene = split[2]
            self.probe_bp = int(split[3])
            self.top_snp = split[4]
            self.top_snp_chr = int(split[5])
            self.top_snp_bp = int(split[6])
            self.allele_1 = split[7]
            self.allele_2 = split[8]
            self.allele_freq = float(split[9])
            self.beta_gwas = float(split[10])
            self.se_gwas = float(split[11])
            self.p_gwas = float(split[12])
            self.beta_eqtl = float(split[13])
            self.se_eqtl = float(split[14])
            self.p_eqtl = float(split[15])
            self.beta_smr = float(split[16])
            self.se_smr = float(split[17])
            self.p_smr = float(split[18])
            try:
                self.p_het = float(split[19])
                self.nsnp_het = float(split[20])
                self.has_het_test = True
            except ValueError:
                self.has_het_test = False
        else: #construct using two MA like objects
            self.beta_gwas = outcome_result.get_beta()
            self.se_gwas = outcome_result.get_se()

            self.beta_eqtl = exposure_result.get_beta()
            self.se_eqtl = exposure_result.get_se()

    def get_beta_outcome(self):
        beta = self.beta_gwas
        return beta

    def get_z_score_outcome(self):
        return self.beta_gwas / self.se_gwas

    def get_beta_exposure(self):
        return self.beta_eqtl

    def get_z_score_exposure(self):
        return self.beta_gwas / self.se_gwas




# probably want to turn this into a method at some time.
def estimate_beta_se_smr(smr_result):
    """
    Determine SMR test effect and standard error.

    :param exposure_data: Exposure estimates which must have the methods get_beta and get_z_score
    :param outcome_data: Outcome estimates which must have the methods get_beta and get_z_score
    :return: tuple of smr beta and smr se
    """

    # from SMR paper.
    t_stat = ((smr_result.get_z_score_exposure() ** 2) * (smr_result.get_z_score_outcome() ** 2)) \
             /  \
             ((smr_result.get_z_score_exposure() ** 2) + (smr_result.get_z_score_outcome() ** 2))

    # this needs a check,
    # checked it with results from Zhy et al. see validate_SMR_estimates.py
    p_value = scipy.stats.chi2.sf(t_stat, 1)
    z_score = scipy.stats.norm.isf(p_value / 2)

    # also from Zhu et al.
    beta_smr = smr_result.get_beta_outcome() / smr_result.get_beta_exposure()
    se_smr = abs(beta_smr / z_score)

    return beta_smr, se_smr, p_value