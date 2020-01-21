

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
Everything was tested on Ubuntu 18.04, we expect MR-link to work on other UNIX bases systems including macOS,
provided tha the following requirements are in place:

Requirements are:
* Python >= 3.6

* pip3 for installing packages

* setuptools for installing ``genome_integration`` you can install setuptools using ``pip3 install setuptools``


* ``gcta64`` and ``plink`` (1.9) are installed and in your PATH as the subprocess module will directly refer to them.

To download ``gcta64`` and ``plink``, please follow the following links.

* `Download GCTA <http://cnsgenomics.com/software/gcta/d>`_

* `Download plink <https://www.cog-genomics.org/plink2/>`_

If you want to check if they are in your path, try: ``which plink`` and ``which gcta64``

Additionally, if you want to test the functionality of the ``genome_integration`` library, you need to have installed the
R language:


Installation
------------
If GCTA and plink are in your ``$PATH``, you can install the genome_integration library with the following command: 

   python3 setup.py build install --user

Which will take about 2 minutes to install. If you want to install the genome_integration library for all users, remove
the ``--user`` option from the command.

Next, you can test the library with:

    python3 setup.py test

Now you have installed the library to run MR-link. Running MR-link is possible using the guide on the
`MR-link page <about_mr_link>`_.


