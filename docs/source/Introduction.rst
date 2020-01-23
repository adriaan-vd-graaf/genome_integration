

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

Everything was tested on Ubuntu 18.04 and macOSX, we assume other distributions of linux work as well.
All analyses can be run on standard hardware.

Requirements are:
    * GCC C and C++ compiler (``gcc``, ``g++``)
    * Python >= 3.6
    * pip3 for installing python packages

And the following python3 packages (which can be installed using pip3)
    * setuptools (``pip3 install setuptools``), which installs setuptools. Used to setup the libary.
    * Python wheel (``pip3 install wheel``)
    * Python PYMC3 (``pip3 install pymc3``) for the _p_ value calibration step

Please make sure that ``gcta64`` and ``plink`` (1.9) are in your PATH as ``genome_integration`` and MR-link will directly
refer to these programs.
* [Download GCTA](http://cnsgenomics.com/software/gcta/)
* [Download plink](https://www.cog-genomics.org/plink2/)

If you want to check if they are in your path, try: ``which plink`` and ``which gcta64``

If you run the tests in the code, you also need to install  R: ``Rscript`` needs to be in your path.

Running an example gene for MR-link will take approximately 10 seconds on a quad core Intel Core i7-7700HQ CPU processor
and require up to 8 GB of RAM.

Running the _p_ value calibration script takes approximately 30 minutes, but is only required once after all the genes are
considered.


Installation
------------
If GCTA and plink are in your ``$PATH``, you can install the genome_integration library with the following command:

   python3 setup.py build install --user

Which will take about 2 minutes to install. If you want to install the genome_integration library for all users, remove
the ``--user`` option from the command.

Next, you can test the library with:

    python3 setup.py test

Now you have installed the library to run MR-link.


