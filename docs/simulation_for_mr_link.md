# Simulation of data that can be used for MR-link

This document describes how genotypes, and phenotypes for causal relationships where simulated in the MR-link manuscript.

We will simulate data in two steps:

1. Simulation of genotypes from non finnish european individuals. (fairly long)

2. Simulation of phenotypes that can be used for causal relationships. (a single python script)

If you already have two genotype file from which you want to simulate phenotypes in two separate cohort 

### Requirements for simulations of genotypes

For the simulation of genotypes, we require the following: 
1. a ~20 Gb of empty disk space and about 20Gb of RAM memory
2. bash shell with awk and wget
3. internet access to download the 1000g chromosome 2, 
4. The HAPgen2 program
5. Python version  >= 3.6, possibly other versions of python3
6. The downloaded genome integration package on your filesystem.
7. plink version 1.9

### Requirements for the simulation of phenotypes

The genome_integration package should be installed.
2 genotype files in plink bed format, 


#Step 1. Simulation of genotypes from 1000 genomes.

Before we start. Please go to the directory relative to where `genome_integration` was downloaded / cloned 
`./mr_link/thousand_genomes_non_finnish_european_subset`.

In this folder, the genotypes for our region of interest will be simulated.

Download chromosome 2 from the 1000 genomes FTP website 
```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
```

Isolate a region 100 - 105 Mb from chromosome 2 from the file that was just downloaded which were downloaded and the header and region was retained using an awk script.

```
awk '{ if( (substr($1, 1, 1) == "#") || (($2 > 100000000) && ($2 < 105000000))) print $0 }' <(zcat ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz) | gzip -c > ALL.chr2.phase3.100.to.105.Mb.1000g.vcf.gz
```

To save disk space, the full file was removed and only the small region was kept. 

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

Interpolation was done using the following python script 
```
python3 interpolate_genetic_map.py
```
    

Haplotypes of biallelic snps were isolated from the individuals which are part of the European superpopulation, 
except for finnish individuals. The list of individuals to include were saved in the file `non_finnish_european_individuals.txt`.
which was made based on the sample table from samples in the `20130606_sample_info.xlsx`.
from `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/`. 

We have made `.haps` and `.legend` files for input into HAPgen2 using the following python script:
```
python3 isolate_variants_make_hapgen_compatible.py
```
 
Using the HAPGEN 2 program genotypes were simulated for 25,000 individuals. 
The `hapgen2` executable should be in this directory, please download it.

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
    

Using plink the files are turned into bed files for later use, in the subsequent simulations. 
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


# STEP 2. simulation of causal relationships

The python script `mr_link/simulate_phenotypes.py` can be used to simulate phenotypes which are causally related in 2 different cohorts, 
which are then by default stored  in the `mr_link/simulated_phenotypes/` directory in a file format that the mr link script 
(`mr_link/MRlink.py`) can understand.

It's necessary to specific many things in the `simulate_phenotypes.py` script:
- 


