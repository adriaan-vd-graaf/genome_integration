# Running MR-link

After [installing](Introduction.md) the `genome_integration` library it is possible to run MR-link.

## Requirements of MR-link

MR-link requires the following for it to run
- Individual level genotypes and phenotypes of the outcome (a complex trait) 
- Summary statistics of the exposure (gene expression)
- Genotype of a sufficiently large reference cohort.

If you do not have these requirments, but want to run MR-link from simulated data, it is possible to simulate that 
[here](simulation_for_mr_link.md). Please set the ensg id option to: `--ensg_id simulated_run`, as simulations
do not represent any gene locations.


## Example code of MR-link

Running MR-link is possible using the following command in the `./mr_link/` directory
```bash
python3 MRlink.py \
    --outcome_bed_file <outcome_bed_file> \
    --outcome_phenotype_file <outcome_phenotype_file> \
    --reference_bed_file <reference_bed_file> \
    --exposure_summary_statitics <summary_statistics_file_of_the_exposure>\
    --ensg_id <ensembl_id_of_the_gene> \
    --temporary_location_prepend <a_location_to_store_temporary_files> \
    --p_val_iv_selection_threshold <p_value_selection_threshold>
    --output_file <file_where_to_output_the_result>
    
```

The files and their formats are described below. 

The other options are described below. 
- `--ensg_id` is the ensembl id used to identify the genomic region from where the 
IVs and the causal relationship are estimated. (set this to `simulated_run` if you're running [simulated data](simulation_for_mr_link.md))
- `--temporary_locatoin_prepend` is a location (directory) where to store temporary files.
- `--p_val_iv_selection_threshold` is the _p_ value used for GCTA-COJO.
- `--output_file` is the file where the result is output (appended) to.


#### Genotype files
Genotypes need to be in the plink bed format. If you do not have genotypes, but want to run MR-link you can simulate 
them [here](simulation_for_mr_link.md)   

#### Phenotype files
Phenotype files need to be a tab separated table with the following columns. 
1. `FID` -- family ID that is matched with the plink genotype column.
2. `IID` -- individual ID that is matched with the plink genotype column.
3. `PHENO` -- Numeric phenotype
First line is the header and will be checked against the following:
```
correct_header= "FID\tIID\tPHENO\n"
```

#### Summary statistic files
Phenotype files need to be a tab separated table with the following columns.
 
1. `CHR` -- chromosome identifyer
2. `POS` -- base pair position
3. `NAME` -- name of the SNP
4. `REF_ALLELE` -- reference allele (usually the major allele)
5. `ALT_ALLELE` -- alternative allele (usually the minor allele)
6. `BETA` -- beta effect size of the marginal estimate
7. `SE` -- standard error of the effect size of the marginal estimate.
10. `MAF` -- allele frequency of the alternative allele (minor allele)
11. `N_OBS` -- number of observations used for the estimate.

first line of the file is the header and will be checked against the following:
```
correct_header = "CHR\tPOS\tNAME\tREF_ALLELE\tEFFECT_ALLELE\tBETA\tSE\tMAF\tN_OBS\n"
```
In the analysis a _p_ value is estimated in the summary statistics using a two tailed T distribution with `N_OBS`-2 
degrees of freedom.
