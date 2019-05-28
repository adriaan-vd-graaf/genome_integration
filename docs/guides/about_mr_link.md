# MR-link usage

After installing the `genome_integration` library it is possible to run MR-link if you follow the 


## Requirements of MR-link

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



#Examples 