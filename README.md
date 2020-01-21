# MR-link and `genome_integration`

## Introduction
This package contains MR-link, a Mendelian randomization (MR) method that efficiently identifies causal relationships 
between gene expression and complex traits, implicitly correcting for unobserved pleiotropy even if there is only a 
single instrumental variable available.

#### Structure of the package: 
`genome_integration` is the library and the folder `./mr_link/` has the programs for an 
implementation of MR-link.

**Extensive documentation of MR-link and the genome integration package can be read in our [readthedocs documentation](https://genome-integration.readthedocs.io/en/latest/)**

Please find more details of the method in our [preprint](https://www.biorxiv.org/content/10.1101/671537v1)

If you want to **contribute**, or have any questions please don't hesitate to file a pull request, or submit an issue in github.

## Requirements

Everything was tested on Ubuntu 18.04 and macOSX, we assume other distributions of linux work as well.
All analyses can be run on non-standard hardware. 

Requirements are: 
- GCC C compiler
- GCC C++ compiler (only required for the p value calibration step)
- Python >= 3.6 
- pip3 for installing python packages

And the following python3 packages (which can be installed using pip3) 
- setuptools (`pip3 install setuptools`), which installs setuptools. Used to setup the libary. 
- Python wheel (`pip3 install wheel`)
- Python PYMC3 (`pip3 install pymc3`) for the p value calibration step

Please make sure that `gcta64` and `plink` (1.9) should be in your PATH as the subprocess module will directly refer to them.
- [Download GCTA](http://cnsgenomics.com/software/gcta/)
- [Download plink](https://www.cog-genomics.org/plink2/)

If you want to check if they are in your path, try: `which plink` and `which gcta64`

If you run the tests in the code, you also need to install  R: `Rscript` needs to be in your path.  

Running an example gene for MR-link will take approximately 10 seconds on a quad core Intel Core i7-7700HQ CPU processor 
and require up to 8 GB of RAM.

Running the p value calibration script takes approximately 30 minutes, but is only required once after all the genes are
considered. 

## Installing:
If all the requirements are met, you can install the genome_integration library with the command in 
the downloaded path
```
python3 setup.py build install --user
```
Installation will take approximately 2 minutes.
 
If you want to install the genome_integration library for all Python users, remove the ``--user`` option 
from the command.

Now that the genome_integration library is installed, we can run MR-link, two examples are described below.
More extensive documentation is available at our [readthedocs documentation](https://genome-integration.readthedocs.io/en/latest/)*

## Running MR-link on two example phenotypes

To run MR-link, please go to the `./mr_link` directory.

The first example will be an analysis where the gene (exposure) has been simulated to have no causal effect effect on the outcome.
 
Running MR-link is possible using the following command:
```bash
# This will run an example of a non-causal gene.
python3 MRlink.py --outcome_bed_file example_genotypes/outcome_cohort \
   --reference_bed example_genotypes/reference_cohort \
   --exposure_summary_statistics example_files/no_causal_effect_exposure_sumstats.txt \
   --outcome_phenotype_file  example_files/no_causal_effect_outcome_pheno.txt   \
   --temporary_location_prepend tmp_loc \
   --p_val_iv_selection_threshold 5e-8 \
   --output_file no_causal_effect_example.txt \
   --ensg_id ENSG00000000000
```

The MR-link output will indicate that the exposure is not causal to the outcome. 
Results can be found in the `no_causal_effect_example.txt` file. 
The stdout of the program will contain the following line with the MR-link result for this example. 
```
Uncalibrated MR-link results: beta: -0.0128, se: 0.10784, p value: 9.06e-01
```
Note that the _p_ value is approximately 0.9. 
The input and output file formats are fully documentated below.

Running the command below will run MR-link using an example data set where the gene was simulated to have a non-null causal effect.
```bash
# This will run an example of a gene with a causal effect.
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
The stdout  will contain the following line with the result for this example. 
```
Uncalibrated MR-link results: beta: 0.4150, se: 0.14734, p value: 4.85e-03
```
The pvalue for the causal effect of the exposure ENSG00000000000 is 4.8e-03. 

A full description of the input and output formats of MR-link is located in our
[readthedocs documentation](https://genome-integration.readthedocs.io/en/latest/about_mr_link.html) 

## MR-link input formats.

MR-link requires the following inputs:
- A summary statistics file of the exposure phenotype. (`--exposure_summary_statistics`) 
- A phenotype file of the output (`-outcome_phenotype_file`)
- Outcome genotypes  (plink bed file format) (`--outcome_bed_file`)
- Reference genotypes (plink bed file format) (`--reference_bed_file`)

There are 2 requirements for the input files: 
- The variants in the outcome genotypes and the reference genotypes need to be the same. 
- There needs to be individual overlap between the outcome genotypes and the outcome phenotype file.

###### Summary statistics input file
Summary statistics are required for the exposure phenotype (`--exposure_summary_statistics`). 

The summary statistic files need to be a tab separated table with the following columns.
1. `CHR` -- chromosome identifier
2. `POS` -- base pair position
3. `NAME` -- name of the SNP
4. `REF_ALLELE` -- reference allele (usually the major allele)
5. `ALT_ALLELE` -- alternative allele (usually the minor allele)
6. `BETA` -- beta effect size of the marginal estimate
7. `SE` -- standard error of the effect size of the marginal estimate.
10. `MAF` -- allele frequency of the alternative allele (alternative allele)
11. `N_OBS` -- number of observations used for the estimate.

First line of the file is the header and will be checked against the following:
```
correct_header = "CHR\tPOS\tNAME\tREF_ALLELE\tEFFECT_ALLELE\tBETA\tSE\tMAF\tN_OBS\n"
```

###### Phenotype input file
MR-link requires the outcome phenotype to be directly measured in the individuals that are present in 
the outcome 

Phenotype files need to be a tab separated table with the following columns.
Phenotypes need to be quantitative and should be corrected for covariates before hand.
 
1. `FID` -- family ID that is matched with the plink genotype column.
2. `IID` -- individual ID that is matched with the plink genotype column.
3. `PHENO` -- Numeric phenotype (please note that MR-link has only been tested for quantitative phenotypes.)

The first line is the header and will be checked against the following:
```
correct_header= "FID\tIID\tPHENO\n"
```
 

### MR-link output format
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


## Calibration of p values

After a first pass of MR-link and if you have at least 100 and preferably at least 1,000 uncalibrated p values for 
different  genes, it is possible to calibrate them using the script located in 
`./mr_link/p_value_calibration.py`.

Running a single p value calibration will take up to 30 minutes, but only has to be performed once at the end of 
an analysis, when MR-link has causal estimates for all your genes.

A full description of the input and output formats for MR-link can be found in [readthedocs](https://genome-integration.readthedocs.io/en/latest/about_mr_link.html)

### An example of p value calibration

After installation of PYMC3 It is possible to run the p value calibration script using the following commands

After a first pass of MR-link and if you have at least 100 and preferably at least 1,000 p values for different  genes,
it is possible to calibrate them using the script located in 
`./mr_link/p_value_calibration.py`.
```bash
#Run this from the ./mr_link/ directory
python3 p_value_calibration.py --input_file example_files/uncalibrated_p_values_example.txt --output_file calibrated_p_values.txt
```
Which will output calibrated _p_ values in the `calibrated_p_values.txt` file, and accept uncalibrated p values from the
`uncalibrated_p_values_example.txt` file. 
###### file formats for p value calibration

The `--input_file` from which p values are calibrated should be a tab separated file with two columns:
1. `exposure_name` the name of the exposure that is run
2. `uncalibrated_p_value` the uncalibrated p value from mr link

The `--output_file` is the same file, but with an extra column appended to it:
1. `exposure_name`
2. `uncalibrated_p_value`
3. `calibrated_p_value` the p value after calibration



## Example Simulation of causal phenotypes.

This section gives an example on how to simulate phenotypes that are causal from a single locus. 
In our manuscript this was used  so that we were able to benchmark MR-link and other widely used MR methods. 
We will simulate phenotypes that are (partly) regulated by the cis locus.

Our manuscript relies on simulations of causal phenotypes that were genetically regulated by one locus.
For the simulation of genotypes, we have used the HAPGEN2 program and we have simulated phenotypes using the script 
described here.

In our implementation three phenotypes are simulated: two (potentially) causal exposures and one outcome phenotype. 
For our purposes, we discard the second exposure and only keep a single exposure and the outcome so that we determine 
how well the tested MR methods are able to deal with pleiotropy. In these analyses, the second exposure is discarded  
and is considered as unobserved pleiotropy.

It is possible to simulate phenotypes (after installation of the `genome_integration` package) using the 
`simulate_phenotypes.py` script in the `./mr_link/` subfolder. If you go enter the subfolder, you can simulate 
phenotypes in the following way:

```bash
# This will simulate phenotypes into the simulated_files/ folder
python3 simulate_phenotypes.py --bed_cohort_1 example_genotypes/exposure_cohort \
                               --bed_cohort_2 example_genotypes/outcome_cohort \
                               --out_prepend simulated_files/example \
                               --exposure_1_causal 0.0 \
                               --exposure_2_causal 0.4 \
                               --n_causal_exposure_1 10 \
                               --n_causal_exposure_2 10 \
                               --overlapping_causal 0 \
                               --phenotypes_to_simulate 10 
```

This will save in 10 simulated phenotypes into the `simulated_files/example*` prepend. 
Per simulation run there are 3 files saved, the summary statistics of the exposure (A file that ends with `exposure_sumstats.txt`) 
and the outcome (A file that ends with `outcome_sumstats.txt`) , and the phenotype of the outcome (A file that ends with `outcome_pheno.txt`).  
resulting into a total of 30 files. The file names contain the parameters of the simulation themselves. 

An explanation of the simulation parameters are:
The parameters for the simulations are:  
- `--bed_cohort_1` The cohort from which the exposure is saved 
- `--bed_cohort_2` The cohort from which the outcome is saved
- `--out_prepend` This is the prepend of the phenotype file
- `--exposure_1_causal` The causal effect of the observed exposure
- `--exposure_2_causal` The causal effect of the unobserved exposure
- `--n_causal_exposure_1` The number of causal variants for the observed exposure
- `--n_causal_exposure_1` The number of causal variants for the observed exposure
- `--overlapping_causal` The number of variants that are the same between the two exposures, otherwise they are in LD 0.25 r^2 < 0.95
-- `--phenotypes_to_simulate` the number of phenotypes that are simulated.

For a description of all phenotype simulation parameters, please see the [documentation](https://genome-integration.readthedocs.io/en/latest/simulation_for_mr_link.html).
We have used these phenotypes for the simulations in our manuscript.  

The output files are of the same formats as the summary statistics files and phenotype files, described above.
These files can be directly used as an input to MR-link.

## Running MR-link on your own data

###### Data requirements before you begin 
If you have individual level genotype and data of your outcome of interest and summary statistics on at least 100 _cis_-regulated exposures 