import numpy as np
import scipy.stats
import subprocess
import os
import warnings
from genome_integration import simulate_mr
from genome_integration import utils
from genome_integration.association import GeneticAssociation

def read_assocs_from_plink_qassoc(assoc_file):

    assocs = {}
    with open(assoc_file, "r") as f:
        f.readline()
        for line in f:
            split = line.split()
            for i in range(len(split)):
                if split[i] == "NA":
                    split[i] = np.nan

            snp_name = split[1]
            tmp_assoc = GeneticAssociation(
                dependent_name="sim_pheno",
                explanatory_name=snp_name,
                n_observations = int(split[3]),
                beta = float(split[4]),
                se = float(split[5]),
                r_squared= float(split[6]),
                chromosome=split[0],
                position=split[3],
                major_allele=None,
                minor_allele=None,
                minor_allele_frequency=None,
                reference_allele=None,
                effect_allele=None
            )
            tmp_assoc.set_p_val(float(split[8]))
            assocs[snp_name] = tmp_assoc

    return assocs

def turn_assocs_into_genetic_associations(assocs, ordered_loci, allele_frequency, sample_sizes):
    #warnings turned off for this, as it's a divide by zero sometimes.
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    z_scores = assocs[:,0] / assocs[:,1]
    warnings.filterwarnings("default")

    p_values = scipy.stats.norm.sf(np.abs(z_scores)) *2

    assocs = {ordered_loci[i].snp_name:
            GeneticAssociation(dependent_name="simulation",
                explanatory_name=ordered_loci[i].snp_name,
                n_observations = sample_sizes[i],
                beta=assocs[i,0],
                se=assocs[i,1],
                r_squared = None,
                chromosome = ordered_loci[i].chromosome,
                position = ordered_loci[i].position,
                major_allele = ordered_loci[i].major_allele,
                minor_allele = ordered_loci[i].minor_allele,
                minor_allele_frequency = allele_frequency[i],
                reference_allele = None,
                effect_allele = None
                )
          for i in range(len(assocs))
          }
    [assocs[ordered_loci[i].snp_name].set_p_val(p_values[i]) for i in range(len(assocs))]
    return assocs


def test_compare_plink_assoc():
    np.random.seed(13289)

    rel_path = '/'.join(('test_resources', 'subset_of_exposure_cohort'))
    if len(__file__.split("/")) > 1:
        plink_loc = "{}/{}".format("/".join(__file__.split("/")[:-1]), rel_path)
    else:
        plink_loc = rel_path

    temp_data = '/'.join(('temp_data', 'plink_file_cojo_test'))

    if len(__file__.split("/")) > 1:
        temp_data = "{}/{}".format("/".join(__file__.split("/")[:-1]), temp_data)

    plinkfile = utils.PlinkFile(plink_loc)
    geno_mat = plinkfile.read_bed_file_into_numpy_array()

    #one causal SNP.
    beta = [0.5, 0.5, -0.4]
    phenotypes = simulate_mr.scale_geno_vec(geno_mat[:,5]) * beta[0]
    phenotypes += simulate_mr.scale_geno_vec(geno_mat[:,7]) * beta[1]
    phenotypes += simulate_mr.scale_geno_vec(geno_mat[:, 100]) * beta[2]
    phenotypes += np.random.normal(size=phenotypes.shape)

    phenotypes -= np.mean(phenotypes)
    phenotypes /= np.std(phenotypes)

    #Write and do the plink association.
    pheno_file = temp_data + "_pheno"
    assoc_file = temp_data + "_assoc"
    with open(pheno_file, "w") as f:
        f.write(f"FID\tIID\tPHENO\n")
        for sample_name, phenotype in zip(plinkfile.fam_data.sample_names, phenotypes):
            sample = plinkfile.fam_data.fam_samples[sample_name]

            f.write(f"{sample.fid}\t{sample.iid}\t{phenotype}\n")

    subprocess.run(["plink",
                    "--bfile", plink_loc,
                    "--assoc", "--allow-no-sex",
                    "--pheno", pheno_file,
                    "--out", assoc_file,
                    ], check=True, stdout=subprocess.DEVNULL, stderr = subprocess.DEVNULL)

    plink_ref_assocs = read_assocs_from_plink_qassoc(assoc_file + ".qassoc")
    own_assocs = np.apply_along_axis(simulate_mr.do_gwas_on_scaled_variants,
                                     axis=0,
                                     arr=geno_mat,
                                     dependent=phenotypes).T

    plink_assocs = np.asarray([ [plink_ref_assocs[x].beta, plink_ref_assocs[x].se]
                                for x in plinkfile.bim_data.snp_names])
    plink_assocs[np.isnan(plink_assocs)] = 0.

    #tolerance is relatively low, plink reports about three sig digits. therefore tolerance is low.
    assert(np.all(np.isclose(own_assocs, plink_assocs, rtol=1e-3, atol = 1e-3)))

    #clean up.
    np.random.seed()
    subprocess.run(["rm", "-f",
                    pheno_file,
                    assoc_file + ".log",
                    assoc_file + ".nosex",
                    assoc_file + ".qassoc",
                    ])


rel_path = '/'.join(('temp_data', ''))
if len(__file__.split("/")) >1:
    test_data_dir = "{}/{}".format("/".join(__file__.split("/")[:-1]), rel_path)
else:
    test_data_dir = rel_path
    
if not os.path.isdir(test_data_dir):
    os.mkdir(test_data_dir)

test_compare_plink_assoc()

