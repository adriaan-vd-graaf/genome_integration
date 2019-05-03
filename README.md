#MR-link and `genome_integration`

##Introduction
This package contains MR-link, a Mendelian randomization (MR) method that efficiently identifies causal relationships 
between gene expression and complex traits, implicitly correcting for unobserved pleiotropy even if there is only as 
single instrumental variable available.

Structure of the package: 
`genome_integration` is the library and `./mr_link/` has the programs for an implementation of MR-link.

Please revisit to read it in our upcoming publication.

####Features of MR-link 
Features of MR-link are:
- Identify robust causal relationships between gene expression and complex quantitative traits. 
- Simulate causal phenotypes from a genetic region, specifying causal effect, causal variants and pleiotropic effects.

####Features of `genome_integration`
`genome_integration` is designed as a library with helper functions and analysis functionality that help the analysis 
 and integration of genomics data.  

Features of genome_integration include:

- Classes to handle single nucleotide polymorphism data 
- Classes to handle genetic associations
- Classes to handle individual samples
- Causal inference using many different popular summary statistics based MR methods. (MR-Egger, IVW, MR-PRESSO, ...)
- Causal inference using MR-link
- Reading genotype files in the `plink` bed format
- Exporting to, Running, and loading often used commands in the `plink` and `gcta64` programs.
- Using the `enrichr` api on a set of genes to enrichment analyses. 


##Requirements
Everything was tested on Ubuntu 18.04, other flavors of linux should also work as well likely also macOS

Requirements are: Python >= 3.6
with the following packages:
numpy, scipy, sklearn, statsmodels, requests (for `enrichr` api), bitarray and plinkio (used as a reference).

Please make sure that `gcta64` and `plink` (1.9) should be in your PATH as the subprocess module will directly refer to them.
- [Download GCTA](http://cnsgenomics.com/software/gcta/#Download)
- [Download plink](https://www.cog-genomics.org/plink2/)

If you want to check if they are in your path, try: `which plink` and `which gcta64`


##installing:
If GCTA and plink are in your `$PATH`, you can install the genome_integration library with the command:
```
python3 setup.py build install test
```

#MR-link example

TODO

#Requirements of MR-link

MR-link requires the following for it to run
- Individual level genotypes and phenotypes of the outcome (a complex trait)
- Summary statistics of the exposure (gene expression)
- Genotype of a sufficiently large reference cohort (Could be the outcome cohort)

####Genotype files
Genotypes need to be in the plink bed format.

####Phenotype files
Phenotype files need to be a tab separated table with the following columns. 
1. `FID` -- family ID that is matched with the plink genotype column.
2. `IID` -- individual ID that is matched with the plink genotype column.
3. `PHENO` -- Numeric phenotype
First line is the header and will be checked against the following:
```
correct_header= "FID\tIID\tPHENO\n"
```

####summary statistic files
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

# Simulations of causality 
TODO

