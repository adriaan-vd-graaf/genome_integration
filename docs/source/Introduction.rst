

Introduction
------------

This package contains MR-link, a Mendelian randomization (MR) method that efficiently identifies causal relationships
between gene expression and complex traits, implicitly correcting for unobserved pleiotropy even if there is only one
instrumental variable available.

The source code of the package can be found at `github <https://github.com/adriaan-vd-graaf/genome_integration>`_

Structure of the package
-------------------------
``genome_integration`` is the library and ``./mr_link/`` has the programs for an implementation of MR-link.

For more information please read our `preprint <https://www.biorxiv.org/content/10.1101/671537v1>`_

If you use MR-link, please cite us.


Requirements
--------------------
Everything was tested on Ubuntu 18.04, we expect MR-link to work on other unix bases systems as well such as macOS,
provided tha the following requirements are in place:

Requirements are:
* Python >= 3.6 with the following packages: numpy, scipy, sklearn, statsmodels, requests (for ``enrichr`` api), bitarray and plinkio (used as a reference).

* ``gcta64`` and ``plink`` (1.9) are installed and in your PATH as the subprocess module will directly refer to them.

To download `gcta64` and ``plink``, please follow the following links.

* `Download GCTA <http://cnsgenomics.com/software/gcta/#Download>`_

* `Download plink <https://www.cog-genomics.org/plink2/>`_

If you want to check if they are in your path, try: ``which plink` and ``which gcta64`


Installation
------------
If GCTA and plink are in your ``$PATH``, you can install the genome_integration library with the command in
the downloaded path

   python3 setup.py build install test

Which will take approximately 2 minutes to install
After which you have installed the library to run MR-link. Running MR-link is possible using the guide on the
`MR-link page <about_mr_link>`_.


