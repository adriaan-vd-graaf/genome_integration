# MR-link and `genome_integration`

## Introduction
This package contains MR-link, a Mendelian randomization (MR) method that efficiently identifies causal relationships 
between gene expression and complex traits, implicitly correcting for unobserved pleiotropy even if there is only a 
single instrumental variable available.

#### Table of contents
 This readme contains the following information
* [Introduction](#Introduction)
* [Requirements](#Requirements)
* [Installation](#Installation)
* [MR-link](#mr-link)
    - [MR-link input formats](#mr-link-input-formats)
    - [MR-link output format](#mr-link-output-format)
* [Calibration of _p_ values](#Calibration-of-p-values)
    - [_p_ value calibration file formats](#_p_-value-calibration-file-formats)
* [Step by step guide to an MR-link analysis](#step-by-step-guide-to-an-mr-link-analysis)
* [Example simulation of causal phenotypes](#Example-simulation-of-causal-phenotypes)

**If something is not clear, or if you want more information, please see our more extensive [readthedocs documentation](https://genome-integration.readthedocs.io/en/latest/)**
 
 
`genome_integration` is the library and the folder `./mr_link/` has the programs for an implementation of MR-link.
Please find more details of the method in our **[preprint](https://www.biorxiv.org/content/10.1101/671537v1)**
If you want to contribute, or have any questions please don't hesitate to file a pull request, or submit an issue in github issue tracker.

## Requirements

Everything was tested on Ubuntu 18.04 and macOSX, we assume other distributions of linux work as well.
All analyses can be run on standard hardware. 

Requirements are: 
- GCC C and C++ compiler (`gcc`, `g++`)
- Python >= 3.6 
- pip3 for installing python packages

And the following python3 packages (which can be installed using pip3) 
- setuptools (`pip3 install setuptools`), which installs setuptools. Used to setup the libary. 
- Python wheel (`pip3 install wheel`)
- Python PYMC3 (`pip3 install pymc3`) for the _p_ value calibration step

Please make sure that `gcta64` and `plink` (1.9) are in your PATH as `genome_integration` and MR-link will directly 
refer to these programs.
- [Download GCTA](http://cnsgenomics.com/software/gcta/)
- [Download plink](https://www.cog-genomics.org/plink2/)

If you want to check if they are in your path, try: `which plink` and `which gcta64`

If you run the tests in the code, you also need to install  R: `Rscript` needs to be in your path.  

Running an example gene for MR-link will take approximately 10 seconds on a quad core Intel Core i7-7700HQ CPU processor 
and require up to 8 GB of RAM.

Running the _p_ value calibration script takes approximately 30 minutes, but is only required once after all the genes are
considered. 

## Installation
If all the requirements are met, you can install the `genome_integration` library with the command in 
the folder where you downloaded / cloned the package
```shell script
python3 setup.py build install --user
```
Installation will take approximately 2 minutes.
 
If you want to install the genome_integration library for all Python users, remove the ``--user`` option 
from the command.

Now that the `genome_integration` library is installed, we can run MR-link, two examples are described below. 

###### Testing the library (optional)
Testing the genome integration library is done with the following command:
```shell script
python3 setup.py test
```
Which should pass all tests. If it doesn't, please submit an issue to the issue tracker. 

## MR-link

To run MR-link, please go to the `./mr_link` directory.

The first example will be an analysis where the gene (exposure) has been simulated to have no causal effect effect on the outcome.
 
Running MR-link is possible using the following command:
```shell script
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
The input and output file formats are fully documented below.

Running the command below will run MR-link using an example data set where the gene was simulated to have a causal effect.
```shell script
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
The _p_ value for the causal effect of the exposure ENSG00000000000 is 4.9e-03. 

A full description of the input and output formats of MR-link is located in our
[readthedocs documentation](https://genome-integration.readthedocs.io/en/latest/about_mr_link.html) 

##### MR-link input formats

MR-link requires the following inputs:
- A summary statistics file of the exposure phenotype. (`--exposure_summary_statistics`) 
- A phenotype file of the output (`--outcome_phenotype_file`)
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
the outcome cohort. 

Phenotype files need to be a tab delimited file with the following columns.
Phenotypes need to be quantitative and should be corrected for covariates before hand.
 
1. `FID` -- family ID that is matched with the plink genotype column.
2. `IID` -- individual ID that is matched with the plink genotype column.
3. `PHENO` -- Numeric phenotype (MR-link has only been tested for quantitative phenotypes)

The first line is the header and will be checked against the following:
```
correct_header= "FID\tIID\tPHENO\n"
```

#### MR-link output format
MR-link outputs a file with two lines, one header and one result line for the gene, fields are tab separated
The header looks like this: `"ensembl_name\tmethod\tbeta\tse\tp_value\tn_ivs\tiv_summary\n"`

The result line contains the following fields:
1. `ensembl_name`: the ensembl gene id that was tested.
2. `method` : The method used for causal inference. Currently this is only `MR-link_uncalibrated`, could be expanded in later versions.
3. `beta`: Is the ridge regression beta (This has been regularized so the direction is informative, but the magnitude may not be)
4. `se`: The standard error of the beta estimate
5. `p_value`: is the _p_ value, based on a T test. This is very conservative and needs to be calibrated when all genes have been run.
6. `n_ivs`:  is the number of ivs identified by GCTA-COJO
7. `iv_summary`: is a comma separated field with information on all the ivs. Per iv this is semicolon separated in the following way: `<snp_name>;<effect_allele>;<beta>;<se>;<minor_allele_frequency>` of the iv. the effect size is conditional on the other ivs.   


## Calibration of _p_ values

After a first pass of MR-link and if you have at least 100 and preferably at least 1,000 uncalibrated _p_ values for 
different  genes, it is possible to calibrate them using the script located in 
`./mr_link/p_value_calibration.py`.

Running a single _p_ value calibration will take up to 30 minutes, but only has to be performed once at the end of 
an analysis, when MR-link has causal estimates for all your genes.

A full description of the input and output formats for MR-link can be found in [readthedocs](https://genome-integration.readthedocs.io/en/latest/about_mr_link.html)

###### An example of _p_ value calibration

After installation of PYMC3 It is possible to run the _p_ value calibration.

After the analysis of the `MRlink.py` script resulting in uncalibrated _p_ values it is possible to calibrate them using the  
`./mr_link/p_value_calibration.py` script.

```shell script
#Run this from the ./mr_link/ directory
python3 p_value_calibration.py --input_file example_files/uncalibrated_p_values_example.txt --output_file calibrated_p_values.txt
```

Which will output calibrated _p_ values in the `calibrated_p_values.txt` file, and accept uncalibrated _p_ values from the
`uncalibrated_p_values_example.txt` file. 

#### _p_ value calibration file formats

The `--input_file` from which _p_ values are calibrated should be a tab separated file with two columns:
1. `exposure_name` the name of the exposure that is run
2. `uncalibrated_p_value` the uncalibrated _p_ value from mr link

The `--output_file` is the same file, but with an extra column appended to it:
1. `exposure_name`
2. `uncalibrated_p_value`
3. `calibrated_p_value` the _p_ value after calibration

###### _p_ value calibration of non-null simulations
If you want to calibrate _p_ values without computing the beta distribution, you specify the alpha and beta parameters 
combined with the `--only_calibrate` option in the following way: 
```shell script
#Run this from the ./mr_link/ directory
python3 p_value_calibration.py \
    --input_file example_files/uncalibrated_p_values_example.txt \
    --output_file calibrated_p_values.txt \
    --only_calibrate \
    --alpha_parameter 3.9703 \
    --beta_parameter 0.6106
```
If you use simulated datasets, it's important to calibrate on the null scenarios, then use these parameters on the
non-null scenarios. 



## Step by step guide to an MR-link analysis

This is the basic workflow for using MR-link.

###### Step 0. Data requirements before you begin 

You require the following (large) data for running a succesful MR-link analysis

- One cohort with individual level genotype and phenotype data
- Gene expression summary statistics of multiple gene expression traits*
- (Recommended) a cohort with reference genotypes of sufficient size (n >= 5,000)   

[*] Only gene expression data is supported, other _cis_-regulated phenotypes can be implemented if there is demand. 
Please open an issue in the issue tracker if you want to run other _cis_-regulated phenotypes.

__Nb.__ The source of these data needs to be from an ethnically matched population (genotypes _and_ summary statistics)

###### Step 1. Format your data
You need to format the following data

- Format your genotypes into the plink bed format.
- Format your outcome phenotypes into the [phenotype format](#Phenotype-input-file) 
- Format your exposure summary statistic files into the [summary statistic format](#Summary-statistics-input-file)

###### Step 2. Run MR-link on all the exposures

After data formatting, it is possible to run MR-link for all your gene expression phenotypes as an exposure.
All options for MR-link are documented [here](https://genome-integration.readthedocs.io/en/latest/about_mr_link.html)

###### Step 3. Run the _p_ value calibration

After all the exposures are run and finished, you can calibrate the _p_ values to their expected distribution.
Every gene saves itself into it's own file, which needs to be formatted into an [input file for _p_ calibration](#file-formats-for-p-value-calibration)
If you have a directory with all the MR-link results, you can make an input file for _p_ value calibration with the following command:

```shell script
filename=your_uncalibrated_p_values.txt
echo -e 'exposure_name\tuncalibrated_p_value' > "$filename"
ls | while read mr_link_file
do 
    tail -n 1 "${mr_link_file}" | awk '{print $1"\t"$5}' 
done >> "$filename"
```
Your input file for _p_ value calibration will be saved in the file that is named `your_uncalibrated_p_values.txt`

The full documentation for _p_ value calibration is documented [here](https://genome-integration.readthedocs.io/en/latest/calibrating_mr_link_p_values.html) 


## Example simulation of causal phenotypes

This section gives an example on how to simulate phenotypes that are causal from a single locus. 
In our manuscript this was used  so that we were able to benchmark MR-link and other widely used MR methods. 
We will simulate phenotypes that are (partly) regulated by the _cis_-locus.

Our manuscript relies on simulations of causal phenotypes that were genetically regulated by one locus.
For the simulation of genotypes, we have used the HAPGEN2 program and we have simulated phenotypes using the script 
described here.

In our implementation three phenotypes are simulated: two (potentially) causal exposures and one outcome phenotype. 
For our purposes, we discard the second exposure and only keep a single exposure and the outcome so that we determine 
how well the tested MR methods are able to deal with pleiotropy.

It is possible to simulate phenotypes (after installation of the `genome_integration` package) using the 
`simulate_phenotypes.py` script in the `./mr_link/` subfolder. If you go enter the subfolder, you can simulate 
phenotypes in the following way:

```shell script
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

This will save in 10 sets of simulated phenotypes into the `simulated_files/example*` prepend. 
Per simulation run there are 3 files saved, the summary statistics of the exposure (a file that ends with `exposure_sumstats.txt`), the outcome (a file that ends with `outcome_sumstats.txt`) 
and the phenotype of the outcome (a file that ends with `outcome_pheno.txt`).  
Resulting in a total of 30 files. The file names contain the parameters of the simulation themselves. 
These files are formatted so that `MRlink.py` accepts them.

The parameters for the simulations are:  
- `--bed_cohort_1` The cohort from which the exposure is saved 
- `--bed_cohort_2` The cohort from which the outcome is saved
- `--out_prepend` This is the prepend of the phenotype file
- `--exposure_1_causal` The causal effect of the observed exposure
- `--exposure_2_causal` The causal effect of the unobserved exposure
- `--n_causal_exposure_1` The number of causal variants for the observed exposure
- `--n_causal_exposure_1` The number of causal variants for the observed exposure
- `--overlapping_causal` The number of variants that are the same between the two exposures, otherwise they are in LD 0.25 r^2 < 0.95
- `--phenotypes_to_simulate` the number of phenotype sets that are simulated.

For a description of all phenotype simulation parameters, please see the [documentation](https://genome-integration.readthedocs.io/en/latest/simulation_for_mr_link.html).
We have used these phenotypes for the simulation section of our manuscript.  

The output files are of the [summary statistic format](#Summary-statistics-input-file) and 
[phenotype format](#Phenotype-input-file). These files can be directly used as an input to MR-link.

After running MR-link on the simulated values, it's important to consider that _p_ value calibration should only be done 
on simulation scenarios with no causal effect. To calibrate the _p_ values of a non-null causal effect, save the 
parameters of the beta distribution and use the `--only_calibrate` option in the `p_value_calibration.py` script for the
scenarios with a real causal effect.

