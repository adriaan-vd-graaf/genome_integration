# Simulation of genotypes and phenotypes for MR-link

This document describes how genotypes, and phenotypes for causal relationships where simulated in the MR-link manuscript.

We will simulate data in two steps:

1. Simulation of genotypes from non finnish european individuals. (fairly long)

2. Simulation of phenotypes that can be used for causal relationships. (shorter)

After running these steps, a user can run MR-link to identify causal relationships.    

If you already have two genotype files from which you want to simulate phenotypes in two separate cohort please go to 
step 2.

#### Requirements for simulations of genotypes

For the simulation of genotypes, we require the following: 
1. *30* Gb of empty disk space and about *20Gb* of RAM memory
2. bash shell with awk and wget
3. internet access to download the 1000g chromosome 2 (1.2Gb), 
4. The HAPGEN2 program
5. Python version  >= 3.6.
6. The downloaded genome integration package on your filesystem.
7. plink version 1.9

#### Requirements for the simulation of phenotypes

The genome_integration package should be installed.
2 genotype files in plink bed format, with the same set of variants.


### Step 1. Simulation of genotypes from 1000 genomes.

Before we start. Please go to the directory relative to where `genome_integration` was downloaded / cloned 
`./mr_link/simulation_of_genotypes/`.

In this folder, the genotypes for our region of interest will be simulated.

Download chromosome 2 from the 1000 genomes FTP website 
```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
```

Isolate a region 100 - 105 Mb from chromosome 2 from the file that was just downloaded which were downloaded and the header and region was retained using an awk script.

```
awk '{ if( (substr($1, 1, 1) == "#") || (($2 > 100000000) && ($2 < 105000000))) print $0 }' <(zcat ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz) | gzip -c > ALL.chr2.phase3.100.to.105.Mb.1000g.vcf.gz
```

To save disk space, the full file will be removed below and only the small region was kept. 

```
rm ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
```

Using another awk script, we make another file with only the snp names and the positions used in the next step for the interpolation.

```
awk '{if(substr($1, 1, 1) != "#") print $1, $2, $3}' <(zcat ALL.chr2.phase3.100.to.105.Mb.1000g.vcf.gz) > genetic_variants_in_region.txt 
```

For these genotypes, the recombination rates were interpolated based on the HAPMAP3 recombination rates 

Download this file using the following command
```
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
```
and unpack only the second chromosome using
```
tar -xvf genetic_map_HapMapII_GRCh37.tar.gz genetic_map_GRCh37_chr2.txt
```

Interpolation is  done using the following python script 
```
python3 interpolate_genetic_map.py
```
    

Haplotypes of biallelic snps were isolated from the individuals which are part of the European superpopulation, 
except for Finnish individuals. The list of individuals to include are saved in the file `non_finnish_european_individuals.txt`.
which was made based on the sample table from samples in the `20130606_sample_info.xlsx`.
from `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/`. 

We have made `.haps` and `.legend` files for input into HAPgen2 using the following python script:
```
python3 isolate_variants_make_hapgen_compatible.py
```
 
Using the HAPGEN 2 program genotypes were simulated for 25,000 individuals. 
The `hapgen2` executable should be in your working directory, please download it if you have not done so already.

```
./hapgen2 -m interpolated_genetic_map.100Mb.to.105Mb.map -l europeans_1000g_chr2_100_to_105_Mb.legend -h europeans_1000g_chr2_100_to_105_Mb.hap -o simulated_geno -int 100000000 105000000 -n 25000 0 -dl 100202391 0 0.1 0.2
```

which are stored in the file prepend
`simulated_geno.*`, where the summary filetype shows the command used. 
Control individuals are only considered in later analysis. chromosome names were set using sed:
```
sed 's/^snp_[0-9]\+ /2 /g' simulated_geno.controls.gen > tmp_output.gen
mv tmp_output.gen simulated_geno.controls.gen
```
    

Using plink the files are converted into bed files for use in the subsequent simulations. 
Filtering for MAF 1%  in the full population, and then divided up into two 5,000 individual cohorts and one 15,000 individual 
which are subsequently not MAF filtered, to ensure variants remain the same across cohorts.

The full dataset is saved into full_simulated.(bed|bim|fam), and is divided up into three cohorts:
exposure, reference, and an outcome cohort (15,000 individuals.).

```
plink  --data simulated_geno.controls --maf 0.01 --make-bed --out all_simulated 
```

Make the three separate plink genotype bed files: `exposure_cohort`, `reference_cohort` and `outcome_cohort`
```
awk '{print $1, $2} all_simulated.fam | sort -R > all_individuals_random.txt'
head -n 5000 all_individuals_random.txt > exposure_cohort_individuals.txt
plink --bfile all_simulated --chr 2 --from-mb 102 --to-mb 103 --keep exposure_cohort_individuals.txt --make-bed --out exposure_cohort   

head -n 10000 | tail -n 5000 all_individuals_random.txt > reference_cohort_individuals.txt
plink --bfile all_simulated --chr 2 --from-mb 102 --to-mb 103 --keep reference_cohort_individuals.txt --make-bed --out reference_cohort

tail -n 15000 all_indivduals_random.txt > outcome_cohort_individuals.txt
plink --bfile all_simulated --chr 2 --from-mb 102 --to-mb 103 --keep outcome_cohort_individuals.txt --make-bed --out outcome_cohort
```
And this finalizes the simulation of genotypes for the cohorts, so that we can simulate causal relationships.


### STEP 2. simulation of causal relationships

The python script `mr_link/simulate_phenotypes.py` can be used to simulate phenotypes which are causally related in 2 different cohorts, 
which are then by default stored  in the `mr_link/simulated_phenotypes/` directory in a file format that the MR-link script 
(`mr_link/MRlink.py`) can understand.


In the following example the exposure plink genotype is named `exposure_cohort.(bed|bim|fam)`, 
the outcome plink genotype files are named  `exposure_cohort.(bed|bim|fam)`

It is possible to simulate phenotypes using the `simulate_phenotypes.py` script, if you run 
```bash
python3 mr_link/simulate_phenotypes.py --help
```
you will read the following:
```
usage: simulate_phenotypes.py [-h] --bed_cohort_1 BED_COHORT_1 --bed_cohort_2
                              BED_COHORT_2 [--out_prepend OUT_PREPEND]
                              [--save_as {numpy,text}]
                              [--exposure_1_causal EXPOSURE_1_CAUSAL]
                              [--exposure_2_causal EXPOSURE_2_CAUSAL]
                              [--n_causal_exposure_1 N_CAUSAL_EXPOSURE_1]
                              [--n_causal_exposure_2 N_CAUSAL_EXPOSURE_2]
                              [--overlapping_causal OVERLAPPING_CAUSAL]
                              [--directional_pleiotropy DIRECTIONAL_PLEIOTROPY]
                              [--phenotypes_to_simulate PHENOTYPES_TO_SIMULATE]

optional arguments:
  -h, --help            show this help message and exit
  --bed_cohort_1 BED_COHORT_1
                        The bed file for cohort 1 (from which the exposure is
                        calculated and saved), alleles need to be harmonized
                        with bed_cohort_2
  --bed_cohort_2 BED_COHORT_2
                        The bed file for cohort 2 (from which the outcome is
                        calculated and saved)alleles need to be harmonized
                        with bed_cohort_1
  --out_prepend OUT_PREPEND
                        The prepend for where to place the file
  --save_as {numpy,text}
                        How to save the results. 'numpy' will be more
                        efficient and saves all data to a single file.'text'
                        will be human readable but save to eight separate
                        files per simulation.The latter can blow up the number
                        of files in the file system.
  --exposure_1_causal EXPOSURE_1_CAUSAL
                        causal effect of exposure 1
  --exposure_2_causal EXPOSURE_2_CAUSAL
                        causal effect of the pleiotropic exposure 2
  --n_causal_exposure_1 N_CAUSAL_EXPOSURE_1
                        number of causal SNPs for exposure 1
  --n_causal_exposure_2 N_CAUSAL_EXPOSURE_2
                        number of causal SNPs for the pleiotropic exposure 2
  --overlapping_causal OVERLAPPING_CAUSAL
                        number of causal SNPs that overlap between exposure 1
                        and exposure 2
  --directional_pleiotropy DIRECTIONAL_PLEIOTROPY
                        If the effects of the second exposure should be
                        directional
  --phenotypes_to_simulate PHENOTYPES_TO_SIMULATE
                        How many phenotypes are saved could take a long time.
                        Results are nonetheless saved after each simulation
```

and you can simulate phenotypes using the default settings by running:
```
python3 simulate_phenotypes.py --bed_cohort_1 exposure_cohort --bed_cohort_2 outcome_cohort
```
Which will save simulated phenotypes in the `./mr_link/simulated_files/` folder.
By default, the exposure phenotype has no causal effect on the outcome, and the unobserved (pleiotropic) exposure does.
You can edit the simulation parameters by changing the options, which are explained in the above help.
 

With these genotypes and  phenotypes simulated and the summary statistics calculated, it's now possible to run MR-link.
Please read the [section](about_mr_link.md) that explains how to run MR-link 
