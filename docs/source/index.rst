.. MR-link and genome_integration documentation master file, created by
   sphinx-quickstart on Tue May 28 18:10:10 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MR-link and genome_integration
==========================================================

.. toctree::
   :maxdepth: 0
   :caption: Contents:

   Introduction

   about_mr_link
   simulation_for_mr_link
   calibrating_mr_link_p_values
   about_genome_integration.md

   modules/variants
   modules/samples
   modules/gene_regions
   modules/association
   modules/causal_inference
   modules/simulate_mr
   modules/resources
   modules/utils


Introduction
------------

This package contains MR-link, a Mendelian randomization (MR) method that efficiently identifies causal relationships
between gene expression and complex traits, implicitly correcting for unobserved pleiotropy even if there is only one
single instrumental variable available.

The source code of the package can be found at `github <https://github.com/adriaan-vd-graaf/genome_integration>`_

``genome_integration`` is the library and ``./mr_link/`` has the programs for an implementation of MR-link.

For more information please read our `preprint <https://www.biorxiv.org/content/10.1101/671537v1>`_



