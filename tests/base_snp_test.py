import genome_integration
import genome_integration.variants as variants
import warnings

# this function tests the add_snp_data class.
import genome_integration.variants.SNP


def test_add_snp_snp_acceptance():
    """
    Tests some snps

    :return:
    """
    # making sure that this does not throw an error when the data is the same.
    snp_a = variants.SNP("rs123", "1", 123, "A", "C", 0.25)
    snp_b = variants.SNP("1:123", "1", 123, "A", "C", 0.25)

    failed = False
    try:
        snp_a.add_snp_data(snp_a)
        snp_a.add_snp_data(snp_b)
    except ValueError:
        failed = True

    assert not failed

    # now making sure that this throws an error when comparing snps.
    snp_c = variants.SNP(("2:123", "2", 123, "A", "C", 0.25))

    failed = False
    try:
        snp_a.add_snp_data(snp_c)
    except ValueError:
        failed = True

    assert failed

    snp_d = variants.SNP("rs126", "1", 123, "A", "C", 0.25)

    failed = False
    try:
        snp_a.add_snp_data(snp_d)
    except ValueError:
        failed = True

    assert failed


# this function tests the add_snp_data class.
def test_update_position():

    snp_a = variants.SNP("rs123", "1", 123, "A", "C", 0.25)
    snp_no_position = variants.SNP("rs123", None, None, "A", "C", 0.25)

    snp_no_position.add_snp_data(snp_a)

    assert snp_no_position.chromosome == snp_a.chromosome
    assert snp_no_position.position == snp_a.position

    snp_with_position = variants.SNP("rs123", "2", 124, "A", "C", 0.25)
    snp_with_position.add_snp_data(snp_a)

    assert snp_with_position.chromosome != snp_a.chromosome
    assert snp_with_position.position != snp_a.position


# this function tests the add_snp_data class.
def test_update_alleles():
    """
    This function is used to ensure the alleles are correct when adding snp data.
    :return:
    """

    """
    Making sure the allele data is added correctly
    """
    snp_a = variants.SNP("rs123", "1", 123, "A", "C", 0.25)
    snp_no_alleles = variants.SNP("rs123", None, None, None, None, 0.25)

    snp_no_alleles.add_snp_data(snp_a)
    assert snp_no_alleles.major_allele == snp_a.major_allele
    assert snp_no_alleles.minor_allele == snp_a.minor_allele


    """
    Making sure the snp data is added when there is only one allele present.
    """
    snp_only_major_allele = variants.SNP("rs123", None, None, "A", None, 0.25)
    snp_only_major_allele.add_snp_data(snp_a)
    assert snp_only_major_allele.minor_allele == snp_a.minor_allele

    snp_only_minor_allele = variants.SNP("rs123", None, None, None, "C", 0.25)
    snp_only_minor_allele.add_snp_data(snp_a)
    assert snp_only_major_allele.minor_allele == snp_a.minor_allele

    """
    Making sure the snp data is added when the alleles are flipped.
    """
    allele_frequency = 0.75
    snp_flipped_alleles = variants.SNP("rs123", "1", 123, "C", "A", allele_frequency)
    snp_flipped_alleles.add_snp_data(snp_a)

    assert snp_flipped_alleles.major_allele == snp_a.major_allele
    assert snp_flipped_alleles.minor_allele == snp_a.minor_allele
    assert snp_flipped_alleles.minor_allele_frequency == (1 - allele_frequency)

    """
    Making sure the snp data is added when there is only a major allele. which is the minor allele in the reference.
    """
    snp_flipped_only_major = variants.SNP("rs123", "1", 123, "C", None, allele_frequency)
    snp_flipped_only_major.add_snp_data(snp_a)

    assert snp_flipped_alleles.major_allele == snp_a.major_allele
    assert snp_flipped_alleles.minor_allele == snp_a.minor_allele
    assert snp_flipped_alleles.minor_allele_frequency == (1 - allele_frequency)


    """
    Making sure the snp data is added when there is only a minor allele. which is the major allele in the reference.
    """
    snp_flipped_only_major = variants.SNP("rs123", "1", 123, None, "A", allele_frequency)
    snp_flipped_only_major.add_snp_data(snp_a)

    assert snp_flipped_alleles.major_allele == snp_a.major_allele
    assert snp_flipped_alleles.minor_allele == snp_a.minor_allele
    assert snp_flipped_alleles.minor_allele_frequency == (1 - allele_frequency)

    """
    Check that it fails when both the alleles are wrong.
    """
    snp_both_wrong_alleles = variants.SNP("rs123", "1", 123, "G", "T", 0.25)

    failed = False
    try:
        snp_both_wrong_alleles.add_snp_data(snp_a)
    except ValueError:
        failed = True
    assert failed

    """
    Check that it fails when only a single allele is wrong.
    """
    snp_both_wrong_alleles = variants.SNP("rs123", "1", 123, "A", "T", 0.25)

    failed = False
    try:
        snp_both_wrong_alleles.add_snp_data(snp_a)
    except ValueError:
        failed = True
    assert failed

    snp_both_wrong_alleles = variants.SNP("rs123", "1", 123, "T", "A", 0.25)

    failed = False
    try:
        snp_both_wrong_alleles.add_snp_data(snp_a)
    except ValueError:
        failed = True
    assert failed


def test_adding_minor_allele_frequency():
    """

    :return:
    """

    """
    Tests if minor allele frequency data is correct.
    """
    snp_a = variants.SNP("rs123", "1", 123, "A", "C", 0.25)

    snp_to_add_maf = variants.SNP("rs123", "1", 123, "A", "C")
    snp_to_add_maf.add_snp_data(snp_a)

    assert snp_to_add_maf.minor_allele_frequency == snp_a.minor_allele_frequency

    """
    Tests if minor allele frequency will be flipped. 
    """
    maf = 0.77
    snp_to_add_maf = variants.SNP("rs123", "1", 123, "C", "A", maf)
    snp_to_add_maf.add_snp_data(snp_a)
    assert snp_to_add_maf.minor_allele_frequency == (1-maf)


    """
    Tests if this will produce a runtime warning as there could be allele frequency differences
    """
    warning_given = False
    snp_to_add_maf = variants.SNP("rs123", "1", 123, "C", "A", 1-maf)

    with warnings.catch_warnings(record=True) as w:

        warnings.simplefilter("always")
        snp_to_add_maf.add_snp_data(snp_a)

        if len(w) == 1:
            warning_given = True

    assert warning_given


def test_overwrite_snp_data():
    """
    Testing if the overwrite works.

    :return:
    """

    """
    Overwrite the snp data, this should not fail, and produce the same values in the class.
    """
    snp_a = variants.SNP("rs123", "1", 123, "A", "C", 0.25)

    snp_to_overwrite = variants.SNP("rs123", "hahaa", 1e6, "HGH", "lala", 1000.0)

    snp_to_overwrite.add_snp_data(snp_a, overwrite=True)

    #make sure everythign is the same.
    assert snp_to_overwrite.snp_name == snp_a.snp_name
    assert snp_to_overwrite.chromosome == snp_a.chromosome
    assert snp_to_overwrite.position == snp_a.position
    assert snp_to_overwrite.major_allele == snp_a.major_allele
    assert snp_to_overwrite.minor_allele == snp_a.minor_allele
    assert snp_to_overwrite.minor_allele_frequency == snp_a.minor_allele_frequency

    assert snp_to_overwrite.has_position_data == snp_a.has_position_data
    assert snp_to_overwrite.has_allele_data == snp_a.has_allele_data
    assert snp_to_overwrite.has_frequency_data == snp_a.has_frequency_data

def test_reading_bim_from_snp_pack():
    rel_path = '/'.join(('test_resources', 'subset_of_exposure_cohort.bim'))

    plink_loc = "{}/{}".format("/".join(__file__.split("/")[:-1]), rel_path)
    variants.BimFile(plink_loc)



test_add_snp_snp_acceptance()
test_update_position()
test_update_alleles()
test_adding_minor_allele_frequency()
test_overwrite_snp_data()
test_reading_bim_from_snp_pack()