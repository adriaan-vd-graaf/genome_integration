#MR-link and `genome_integration`

##Introduction
This package contains MR-link, a Mendelian randomization (MR) method that efficiently identifies causal relationships 
between gene expression and complex traits, implicitly correcting for unobserved pleiotropy even if there is only as 
single instrumental variable available.

Structure of the package: 
`genome_integration` is the library and `./mr_link/` has the programs for an implementation of MR-link.

Please revisit to read it in our upcoming publication.


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
If GCTA and plink are in your `$PATH`, you can install the genome_integration library with the command in 
the downloaded path
```
python3 setup.py build install test
```
After which you can run MR-link, details of which are on the MR-link page.