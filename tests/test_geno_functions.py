from genome_integration import simulate_mr
import numpy as np
import scipy.stats

def test_allele_frequency():

    n_size = 10000
    sim_frq = 0.3
    #this is the plinkio encoding
    np.random.seed(24994993)
    genotype_vector = 2 - np.random.binomial(2, sim_frq, n_size)
    geno_frq = simulate_mr.geno_frq(genotype_vector)

    assert(np.isclose(geno_frq, 0.3012))

    #introduce some 3 values (plinkio encodes the missing values as 3)
    missing_num = 100
    genotype_vector[np.random.choice(n_size, missing_num, replace=False)] = 3.0
    geno_frq = simulate_mr.geno_frq(genotype_vector)
    assert (np.isclose(geno_frq, 0.3019696))


def test_scaling():
    n_size = 10000
    sim_frq = 0.01
    # this is the plinkio encoding
    genotype_vector = 2 - np.random.binomial(2, sim_frq, n_size)

    assert(np.isclose(np.mean(simulate_mr.scale_geno_vec(genotype_vector)), 0.0))

    # introduce some 3 values (plinkio encodes the missing values as 3)
    missing_num = 9000
    genotype_vector[np.random.choice(n_size, missing_num, replace=False)] = 3.0
    scaled_geno = simulate_mr.scale_geno_vec(genotype_vector, remove_three_encoding=True)
    assert(len(np.unique(scaled_geno[scaled_geno != 0.0])) <= 3)



test_allele_frequency()
test_scaling()
