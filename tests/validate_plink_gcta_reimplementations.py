import numpy as np
import scipy.stats
import subprocess
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


def test_compare_gcta_cojo_with_reference_no_missing():

    np.random.seed(1313)

    rel_path = '/'.join(('test_resources', 'subset_of_exposure_cohort'))
    plink_loc = "{}/{}".format("/".join(__file__.split("/")[:-1]), rel_path)

    temp_data = '/'.join(('temp_data', 'gcta_cojo_test'))
    temp_data = "{}/{}".format("/".join(__file__.split("/")[:-1]), temp_data)

    plinkfile = utils.PlinkFile(plink_loc)
    geno_mat = plinkfile.read_bed_file_into_numpy_array()
    print(geno_mat.shape)
    #one causal SNP.
    # one causal SNP.
    beta = [2, 2, -4]
    phenotypes = simulate_mr.scale_geno_vec(geno_mat[:, 25]) * beta[0]
    phenotypes += simulate_mr.scale_geno_vec(geno_mat[:, 7]) * beta[1]
    phenotypes += simulate_mr.scale_geno_vec(geno_mat[:, 100]) * beta[2]
    phenotypes += np.random.normal(size=phenotypes.shape)

    phenotypes -= np.mean(phenotypes)
    phenotypes /= np.std(phenotypes)

    assocs = np.apply_along_axis(simulate_mr.do_gwas_on_scaled_variants, axis=0, arr=geno_mat, dependent=phenotypes).T


    minor_allele_frequency = np.apply_along_axis(simulate_mr.geno_frq ,axis=0, arr=geno_mat).T
    ordered_loci = np.asarray([plinkfile.bim_data.bim_results[x] for x in plinkfile.bim_data.snp_names])

    sample_sizes = np.sum(geno_mat != 3, axis=0)
    genetic_assocs = turn_assocs_into_genetic_associations(assocs,
                                                           ordered_loci,
                                                           minor_allele_frequency,
                                                           sample_sizes)



    gcta_assocs = utils.do_gcta_cojo_joint_on_genetic_associations(genetic_assocs,
                                                        plink_loc,
                                                        temp_data,
                                                        p_val_thresh=0.05,
                                                        maf= 0.0, _keep_ma_files=True)

    #update the associations with the correct precision, GCTA also uses.
    with open(temp_data + "_temp.ma", "r") as f:
        f.readline()
        for line in f:
            split = line.split()
            snp_name = split[0]
            frq = float(split[3])
            beta = float(split[4])
            se = float(split[5])
            genetic_assocs[snp_name].beta, \
                genetic_assocs[snp_name].se, \
                genetic_assocs[snp_name].minor_allele_frequency = (beta, se, frq)

    subprocess.run(["rm", temp_data + "_temp.ma"], check=True)

    selected_variants = sorted(gcta_assocs.ma_results.keys())

    gcta_effects = np.asarray([[gcta_assocs.ma_results[x].beta, gcta_assocs.ma_results[x].se]
                               for x in selected_variants])

    indices = np.asarray([plinkfile.bim_data.snp_names.index(x) for x in selected_variants], dtype=int)

    assocs_of_interest = {genetic_assocs[x].explanatory_name: genetic_assocs[x]
                    for x, i in zip(plinkfile.bim_data.snp_names, range(len(ordered_loci))) if i in indices}


    gcta_results_own = simulate_mr.do_gcta_cojo_conditional(
        np.apply_along_axis(simulate_mr.scale_geno_vec, 0, geno_mat[:,indices]),
        assocs_of_interest,
        np.arange(3),
        ordered_loci[indices]
        )

    # Effects deviate slightly due to floating point errors and inverting.
    # last digits deviate slightly.
    print(np.concatenate((gcta_results_own, gcta_effects), axis=1))
    assert(np.all(np.isclose(gcta_results_own, gcta_effects, rtol=1e-4, atol=1e-3)))





def test_compare_plink_assoc():
    np.random.seed(13289)

    rel_path = '/'.join(('test_resources', 'subset_of_exposure_cohort'))
    plink_loc = "{}/{}".format("/".join(__file__.split("/")[:-1]), rel_path)

    temp_data = '/'.join(('temp_data', 'plink_file_cojo_test'))
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



test_compare_plink_assoc()
test_compare_gcta_cojo_with_reference_no_missing()
