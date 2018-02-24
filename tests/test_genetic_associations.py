from genome_integration import association
from genome_integration import variants

def test_adding_snp_to_genetic_association():

    beta=0.1
    se = 0.05
    snp_a = variants.BaseSNP("rs123", "1", 123, "A", "C", 0.25)
    assoc_a = association.GeneticAssociation(
        "CeD",
        "rs123",
        200,
        beta,
        se,
        r_squared=None,
        chromosome=None,
        position=None,
        major_allele=None,
        minor_allele=None,
        minor_allele_frequency=None,
        reference_allele=None,
        effect_allele=None
    )

    assoc_a.add_snp_data(snp_a)

    assert assoc_a.snp_name == snp_a.snp_name
    assert assoc_a.chromosome == snp_a.chromosome
    assert assoc_a.position == snp_a.position
    assert assoc_a.major_allele == snp_a.major_allele
    assert assoc_a.minor_allele == snp_a.minor_allele
    assert assoc_a.minor_allele_frequency == snp_a.minor_allele_frequency

    assert assoc_a.has_position_data == snp_a.has_position_data
    assert assoc_a.has_allele_data == snp_a.has_allele_data
    assert assoc_a.has_frequency_data == snp_a.has_frequency_data


    assert beta == assoc_a.beta
    assert beta / se == assoc_a.z_score

    assoc_b = association.GeneticAssociation(
        "CeD",
        "rs123",
        200,
        beta,
        se,
        r_squared=None,
        chromosome=None,
        position=None,
        major_allele="C",
        minor_allele="A",
        minor_allele_frequency=None,
        reference_allele=None,
        effect_allele=None
    )

    assoc_b.add_snp_data(snp_a)

    assert beta == (-1 * assoc_b.beta)
    assert beta / se == (-1 * assoc_b.z_score)

test_adding_snp_to_genetic_association()