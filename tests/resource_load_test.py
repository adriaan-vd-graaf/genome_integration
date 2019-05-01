from genome_integration.resources import read_gene_information

def test_read_gene_info():
    ensembl_info = read_gene_information()
    ensembl_info.get_sorted_genes()
    assert ensembl_info.str_to_gene("ENSG00000115607") == "IL18RAP"