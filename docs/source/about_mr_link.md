# Running MR-link

After [installing](Introduction.md) the `genome_integration` library it is possible to run MR-link.


## Requirements of MR-link

MR-link requires the following for it to run
- Individual level genotypes and phenotypes of the outcome (a complex trait) 
- Summary statistics multiple exposures (gene expression)
- Genotype of a sufficiently large reference cohort (>= 5000 individuals).

Running an example gene for MR-link will take approximately 10 seconds on a quad core Intel Core i7-7700HQ CPU processor 
and require up to 8 Gb of RAM.

If you want to simulate your own genotypes and phenotypes, it is possible to simulate that 
[here](simulation_for_mr_link.md). Please set the ensg_id option to: `--ensg_id simulated_run`, as simulations
do not represent any gene locations.


## Quick examples to run MR-link

To run MR-link, please go to the `./mr_link` directory.
Running MR-link is possible using the following command:
```bash
# This will run a single example gene that shows no causal effect
python3 MRlink.py --outcome_bed_file example_genotypes/outcome_cohort \
   --reference_bed example_genotypes/reference_cohort \
   --exposure_summary_statistics example_files/no_causal_effect_exposure_sumstats.txt \
   --outcome_phenotype_file  example_files/no_causal_effect_outcome_pheno.txt   \
   --temporary_location_prepend tmp_loc \
   --p_val_iv_selection_threshold 5e-8 \
   --output_file no_causal_effect_example.txt \
   --ensg_id ENSG00000000000
```

Which will show MR-link output for a simulation scenario where the exposure is not causal to the outcome. 
Results will be in the `no_causal_effect_example.txt` file. 
Output will contain the following  line with the MR-link result for this example. 
Note that the p value is approximately 0.9 
```
Uncalibrated MR-link results: beta: -0.0128, se: 0.10784, p value: 9.06e-01
```

Running the command below will run MR-link with a causal effect.

```bash
# This will a command with a causal effect
python3 MRlink.py --outcome_bed_file example_genotypes/outcome_cohort \
   --reference_bed example_genotypes/reference_cohort \
   --exposure_summary_statistics example_files/yes_causal_effect_exposure_sumstats.txt \
   --outcome_phenotype_file  example_files/yes_causal_effect_outcome_pheno.txt   \
   --temporary_location_prepend tmp_loc \
   --p_val_iv_selection_threshold 5e-8 \
   --output_file yes_causal_effect_example.txt \
   --ensg_id ENSG00000000000
```
Results will be in the `yes_causal_effect_example.txt` file.
The standard out will contain the  following line with the result for this example. 
```
Uncalibrated MR-link results: beta: 0.4150, se: 0.14734, p value: 4.85e-03
```

Below, we discuss the `MRlink.py` script in more details.

The _p_ values of the MR-link results are very conservative, therefore we have calibrated our p values using a beta  
distribution. Instructions to calibrate p values can be found [here](calibrating_mr_link_p_values.md).

## MR-link specifications.

Running MR-link is possible using the following command in the `./mr_link/` directory
```bash
python3 MRlink.py \
    --outcome_bed_file <outcome_bed_file> \
    --outcome_phenotype_file <outcome_phenotype_file> \
    --reference_bed_file <reference_bed_file> \
    --exposure_summary_statitics <summary_statistics_file_of_the_exposure>\
    --ensg_id <ensembl_id_of_the_gene> \
    --temporary_location_prepend <a_location_to_store_temporary_files> \
    --p_val_iv_selection_threshold <p_value_selection_threshold> \ 
    --output_file <file_where_to_output_the_result>
    
```

The files and their formats are described below. 

Other options for MR-link are:
- `--ensg_id` is the ensembl id used to identify the genomic region from where the 
IVs and the causal relationship are estimated. (set this to `simulated_run` if you're running the examples above or 
[simulated data](simulation_for_mr_link.md))
- `--temporary_locatoin_prepend` is a location (directory) where to store temporary files.
- `--p_val_iv_selection_threshold` is the _p_ value used for GCTA-COJO.
- `--output_file` is the file where the result is output (appended) to.

#### MRlink.py output file.

MR-link outputs a file with two lines, one header and one result line for the gene, results are tab separated
The header looks like this: `"ensembl_name\tmethod\tbeta\tse\tp_value\tn_ivs\tiv_summary\n"`

Then the result line contains the following fields:

1. `ensembl_name`: the ensembl gene id that was tested.
2. `method` : The method used for causal inference. Currently this is only `MR-link_uncalibrated`, could be expanded in later versions.
3. `beta`: Is the ridge regression beta (this has been regularized so the direction is more important than the magnitude)
4. `se`: The standard error of the beta estimate
5. `p_value`: is the _p_ value, based on a T test. This is very conservative and needs to be calibrated when all genes have been run.
6. `n_ivs`:  is the number of ivs identified by GCTA-COJO
7. `iv_summary`: is a comma separated field with information on all the ivs. Per iv this is semicolon separated in the following way: `<snp_name>;<effect_allele>;<beta>;<se>;<minor_allele_frequency>` of the iv. the effect size is conditional on the other ivs.   

#### Genotype files

Genotypes need to be in the plink bed format, and should contain exactly the same variants in both the reference and the outcome cohort. 
If you do not have genotypes, but want to run MR-link you can simulate them [here](simulation_for_mr_link.md)   

#### Phenotype files
Phenotype files need to be a tab separated table with the following columns.
Phenotypes need to be quantitative and should be corrected for covariates before hand.
 
1. `FID` -- family ID that is matched with the plink genotype column.
2. `IID` -- individual ID that is matched with the plink genotype column.
3. `PHENO` -- Numeric phenotype (please note that MR-link has only been tested for quantitative phenotypes.)

The first line is the header and will be checked against the following:
```
correct_header= "FID\tIID\tPHENO\n"
```

#### Summary statistic file
Phenotype files need to be a tab separated table with the following columns.
 
1. `CHR` -- chromosome identifier
2. `POS` -- base pair position
3. `NAME` -- name of the SNP
4. `REF_ALLELE` -- reference allele (usually the major allele)
5. `ALT_ALLELE` -- alternative allele (usually the minor allele)
6. `BETA` -- beta effect size of the marginal estimate
7. `SE` -- standard error of the effect size of the marginal estimate.
10. `MAF` -- allele frequency of the alternative allele (alternative allele)
11. `N_OBS` -- number of observations used for the estimate.

first line of the file is the header and will be checked against the following:
```
correct_header = "CHR\tPOS\tNAME\tREF_ALLELE\tEFFECT_ALLELE\tBETA\tSE\tMAF\tN_OBS\n"
```
