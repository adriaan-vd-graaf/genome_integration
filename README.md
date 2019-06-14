# MR-link and `genome_integration`

## Introduction
This package contains MR-link, a Mendelian randomization (MR) method that efficiently identifies causal relationships 
between gene expression and complex traits, implicitly correcting for unobserved pleiotropy even if there is only as 
single instrumental variable available.

Structure of the package: 
`genome_integration` is the library and `./mr_link/` has the programs for an implementation of MR-link.

*More extensive documentation can be read at [readthedocs](https://genome-integration.readthedocs.io/en/latest/)*

Please revisit to read it in our upcoming publication.

## Requirements
Everything was tested on Ubuntu 18.04, other flavors of linux should also work as well likely also macOS

Requirements are: Python >= 3.6

Please make sure that `gcta64` and `plink` (1.9) should be in your PATH as the subprocess module will directly refer to them.
- [Download GCTA](http://cnsgenomics.com/software/gcta/#Download)
- [Download plink](https://www.cog-genomics.org/plink2/)

If you want to check if they are in your path, try: `which plink` and `which gcta64`

Running an example gene for MR-link will take approximately 10 seconds on a quad core Intel Core i7-7700HQ CPU processor 
and require up to 8 GB of RAM.  Running permutations take approximately 3 minutes per gene, but is not necessary 
for non-significant genes.

## Installing:
If GCTA and plink are in your `$PATH`, you can install the genome_integration library with the command in 
the downloaded path
```
python3 setup.py build install test
```
After which you can run MR-link, details of which are below, or more extensively described in 
[readthedocs](https://genome-integration.readthedocs.io/en/latest/)*

## Running MR-link on example data

To run MR-link, please go to the `./mr_link` directory.
Running MR-link is possible using the following command:

```bash
python3 MRlink.py --outcome_bed_file example_genotypes/outcome_cohort \
   --reference_bed example_genotypes/reference_cohort \
   --exposure_summary_statistics example_files/no_causal_effect_exposure_sumstats.txt \
   --outcome_phenotype_file  example_files/no_causal_effect_outcome_pheno.txt   \
   --temporary_location_prepend tmp_loc \
   --p_val_iv_selection_threshold 5e-8 \
   --output_file no_causal_effect_example.txt \
   --ensg_id simulated_run \
   --permute False
```

Which will show MR-link output for a simulation scenario where the exposure is not causal to the outcome. 
Results will be in the `no_causal_effect_example.txt` file.

Running the command below will run MR-link with a causal effect.

```bash
python3 MRlink.py --outcome_bed_file example_genotypes/outcome_cohort \
   --reference_bed example_genotypes/reference_cohort \
   --exposure_summary_statistics example_files/yes_causal_effect_exposure_sumstats.txt \
   --outcome_phenotype_file  example_files/yes_causal_effect_outcome_pheno.txt   \
   --temporary_location_prepend tmp_loc \
   --p_val_iv_selection_threshold 5e-8 \
   --output_file yes_causal_effect_example.txt \
   --ensg_id simulated_run \
   --permute False
```
Results will be in the `yes_causal_effect_example.txt` file.


## Calibration of p values

After a first pass of MR-link and if you have at least 100 and preferably 1,000 uncalibrated p values for different 
genes, it is possible to calibrate them using the script located in `./mr_link/p_value_calibration.py` script.

for this you require the `PYMC3` package. You can install this package using
``` bash
pip3 install pymc3
```
Running a single p value calibration will take up to 30 minutes, but only has to be performed once at the end of 
an analysis, when all the genes are run.


### p value calibration example

To run p value calibration you need to install pymc3:

```bash
pip3 install pymc3
```

After installation of PYMC3 It is possible to run the p value calibration script using the following commands

```bash
#Run this from the ./mr_link/ directory
python3 p_value_calibration.py --input_file example_files/uncalibrated_p_values_example.txt --output_file calibrated_p_values.txt
```
Which will output calibrated p values in the `calibrated_p_values.txt` file, and accept uncalibrated p values from the
`uncalibrated_p_values_example.txt` file.