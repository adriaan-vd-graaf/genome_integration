# The genome_integration package

This is a python library used to integrate specific genetic information.
This library comes with the MR-link software, which can be used for causal inference of gene expression traits to 
complex traits


## Features of genome_integration
`genome_integration` is designed as a library with helper functions and analysis functionality that aids in the analysis 
 and integration of genomics data.  

Features of genome_integration include:

- Classes to handle single nucleotide polymorphism data 
- Classes to handle genetic associations
- Classes to handle individual samples
- Causal inference using many different popular summary statistics based MR methods. 
     
     1. IVW,
     2. IVW combined with the q test. 
     3. MR-Egger, 
     4. LDA-MR-Egger,
     4. MR-PRESSO
     
- Causal inference using MR-link
- Reading genotype files in the `plink` bed format
- Exporting to, running, and loading often used commands in the `plink` and `gcta64` programs.
- Using the `enrichr` web api on a set of genes to enrichment analyses. 


## Examples of usage for genome_integration

##Harmonization of associations 
Harmonization of associations is important if you have multiple traits or cohorts with associations, and you want to compare them between one another. Alleles need to be identical to have identical effect alleles summary statistics for comparison. In an MR context harmonization is necessary so that we are we are i) 
comparing the effect size for the same allele, 2) check if inconsistencies in alleles are present for the same data 
set. 

It is possible to harmonize alleles using the GeneticAssociation classes. This will flip betas (multiply by -1), and take the inverse  of the allele frequency if effect alleles are not consitent.

```python
from genome_integration import association

exposure_association = association.GeneticAssociation(
                                dependent_name="example_exposure",
                                explanatory_name="example_variant_name",
                                n_observations = 5000,
                                beta=0.5,
                                se=0.1,
                                r_squared = None,
                                chromosome = 6,
                                position = 1000000,
                                major_allele = "A",
                                minor_allele = "C",
                                minor_allele_frequency = 0.1
                                )



outcome_association = association.GeneticAssociation(
                                dependent_name="example_outcome",
                                explanatory_name="example_variant_name",
                                n_observations = 15000,
                                beta=0.2,
                                se=0.05,
                                r_squared = None,
                                chromosome = 6,
                                position = 1000000,
                                major_allele = "C",
                                minor_allele = "A",
                                minor_allele_frequency = 0.89
                                )

#harmonize to the exposure
outcome_association.add_snp_data(exposure_association)

#major allele is now A for the outcome association.
print(outcome_association.major_allele)
#allele frequency is now 1-(old)
print(outcome_association.minor_allele_frequency)
#beta is flipped.
print(outcome_association.beta)
```

### Mendelian randomization

`genome_integration` Is able to do MR analyses in the following way: 

```python
from genome_integration import causal_inference
#eQTL betas, exposure tuples are beta and SE of the estimates.
exposure_tuples = [ (0.5, 0.1),
                    (0.3, 0.1),
                    (0.2, 0.1)
                    ]
#outcome tuples same format
outcome_tuples = [ (0.2, 0.05),
                   (0.12, 0.07),
                   (0.1, 0.04)
                    ]

mr_class = causal_inference.MendelianRandomization()
for exposure_tuple, outcome_tuple in zip(exposure_tuples, outcome_tuples):
    mr_class.do_and_add_smr_estimation(exposure_tuple, outcome_tuple, None, None, None)

#ivw
mr_class.do_ivw_estimation()


#MR-Egger
mr_class.do_egger_regression()
```
