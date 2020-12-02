import os
import sklearn
import subprocess
import bgen_reader
import copy
import warnings

from genome_integration.variants import SNP
from genome_integration.utils import PlinkFile
from .. samples import Sample
from .gcta_ma_utils import *


class BgenReader_using_library:
    """
    This class is used to read in bgen files, using the bgen_reader library.

    Currently mostly untested.
    Your mileage may vary.

    """
    def __init__(self, bgen_location, sample_location=None):

        """
        This will initialize an opened bgen file.
        It will not read any genotypes, but it will read individuals and variant properties.


        :param bgen_location:
        :param sample_location:
        """

        self.bgen_location = bgen_location
        self.sample_location = sample_location
        self.performed_variant_selection = False
        self.geno_expectations = None

        if not os.path.exists(self.bgen_location):
            raise ValueError(f"No file present in {self.bgen_location=}")
        if not os.path.exists(self.sample_location):
            raise ValueError(f"No file present in {self.sample_location=}")

        self.bgen_handle = bgen_reader.open_bgen(self.bgen_location, self.sample_location, verbose=False)

        #Only bi-allelic variants are supported
        self.bi_allelic_variants = self.bgen_handle.nalleles == 2

        if np.sum(self.bi_allelic_variants) != len(self.bgen_handle.nalleles):
            print("Warning: bgen file contains non-biallelic variants, skipping them as they are not implemented")

        self.variant_names = self.bgen_handle.rsids[self.bi_allelic_variants]
        allele_splits = [allele.split(',') for allele in self.bgen_handle.allele_ids[self.bi_allelic_variants]]

        #turning it into the SNP-dict format as it eases adding information and harmonizing alleles.

        self.variant_dict = {snp_name : SNP(snp_name=snp_name,
                                 chromosome=chromosome,
                                 position=position,
                                 major_allele=allele_splits[1],
                                 minor_allele=allele_splits[0],
                                 minor_allele_frequency=None) for
                             snp_name, chromosome, position, allele_splits in zip(self.variant_names,
                                                                                  self.bgen_handle.chromosomes[
                                                                                      self.bi_allelic_variants],
                                                                                  self.bgen_handle.positions[
                                                                                      self.bi_allelic_variants],
                                                                                  allele_splits
                                                                                )
                             }

        self.variant_dict_pos_name = {f'{chromosome}:{position}': SNP(snp_name=snp_name,
                                           chromosome=chromosome,
                                           position=position,
                                           major_allele=allele_splits[1],
                                           minor_allele=allele_splits[0],
                                           minor_allele_frequency=None) for
                                    snp_name, chromosome, position, allele_splits in zip( self.variant_names,
                                                                              self.bgen_handle.chromosomes[
                                                                                  self.bi_allelic_variants],
                                                                              self.bgen_handle.positions[
                                                                                  self.bi_allelic_variants],
                                                                              allele_splits
                                                                              )
                                     }

        #sample names
        if self.sample_location is None:
            self.sample_names = [f'{x}~__~{x}' for x in self.bgen_handle.samples]
        else:
            self.sample_names = [f'{x.split()[0]}~__~{x.split()[1]}' for x in self.bgen_handle.samples]


    def isolate_region_and_variants_return_genotype_expectation(self, region_to_keep=None, variant_names_to_keep=None, warn=True):
        """

        :param region_to_keep:
        :param variant_names_to_keep:
        :return:
        """

        if self.performed_variant_selection:
            raise NotImplementedError("More than one variant selection has not been implemented")

        if region_to_keep is None and variant_names_to_keep is None and warn:
            print(
                "Warning, no filters specified. If the genotype matrix is large, this can take a lot of memory and time.")
            print(
                f"    Will now try and output a matrix of size {self.bgen_handle.shape}")


        rolling_filter = self.bi_allelic_variants
        if region_to_keep is not None:
            chromosome_filter = region_to_keep.chromosome == self.bgen_handle.chromosomes
            position_filter = np.logical_and(self.bgen_handle.positions >= region_to_keep.start, self.bgen_handle.positions <= region_to_keep.end)
            region_filter = np.logical_and(chromosome_filter, position_filter)
            rolling_filter =np.logical_and(rolling_filter, region_filter)
        if variant_names_to_keep is not None:
            variant_filter = np.zeros(rolling_filter.shape, dtype=bool)
            indices = np.asarray(
                [x for x in range(len(self.variant_names)) if self.variant_names[x] in variant_names_to_keep], dtype=int)
            variant_filter[indices] = True
            rolling_filter = np.logical_and(rolling_filter, variant_filter)

        self.geno_expectations = self.bgen_handle.allele_expectation(rolling_filter)[:, :, 0]


        #remove unused variants from the object. (they will remain in the bgen handle.)

        self.variant_names = self.bgen_handle.rsids[rolling_filter]
        allele_splits = [allele.split(',') for allele in self.bgen_handle.allele_ids]

        # turning it into the SNP-dict format as it eases adding information and harmonizing alleles.

        self.variant_dict = {snp_name: SNP(snp_name=snp_name,
                                           chromosome=chromosome,
                                           position=position,
                                           major_allele=allele_splits[1],
                                           minor_allele=allele_splits[0],
                                           minor_allele_frequency=None) for
                             i, snp_name, chromosome, position, allele_splits in
                             enumerate(zip(self.bgen_handle.rsids,
                                           self.bgen_handle.chromosomes[
                                              self.bi_allelic_variants],
                                           self.bgen_handle.positions[
                                              self.bi_allelic_variants],
                                           allele_splits
                                           )) if rolling_filter[i]
                             }

        self.variant_dict_pos_name = {f'{chromosome}:{position}': SNP(snp_name=snp_name,
                                                                      chromosome=chromosome,
                                                                      position=position,
                                                                      major_allele=allele_splits[1],
                                                                      minor_allele=allele_splits[0],
                                                                      minor_allele_frequency=None) for
                                      i, snp_name, chromosome, position, allele_splits in
                                      enumerate(zip(self.bgen_handle.rsids,
                                                    self.bgen_handle.chromosomes[
                                                        self.bi_allelic_variants],
                                                    self.bgen_handle.positions[
                                                        self.bi_allelic_variants],
                                                    allele_splits
                                                    )) if rolling_filter[i]
                                      }

        self.performed_variant_selection = True

        return self.geno_expectations


    def harmonize_genotypes(self, other_bgen_file, other_genotypes ):
        """
        Harmonizes the other plink file to the alleles of self.
        requires that the variants are the same between files.
        requires that the variants have the same alleles

        WARNING: Will mean that the other plink file will have flipped major and minor alleles.

        :param other_plink_file: Plink File to harmonize
        :return: PlinkFile
        """

        if self.geno_expectations is None and other_bgen_file.geno_expectations is None:
            raise ValueError("Need to perform genotype isolation first on both objects")

        if self.geno_expectations is None:
            raise ValueError("Need to perform genotype isolation first on this object")

        if other_bgen_file.geno_expectations is None:
            raise ValueError("Need to perform genotype isolation on the other object before continueing")

        if set(self.variant_names) != set(other_bgen_file.variant_names):
            raise ValueError("Both genotype files need to have exactly the same set of variants.")

        if self.variant_names != other_bgen_file.variant_names:
            #This is a runtime check to ensure that the variants are ordered the same.
            #One assumption is that plink orders the variants bades on position,
            #so this should not be reached, otherwise ordering functionality needs to be implemented.
            raise ValueError("Ordering of the variants needs to be exactly the same.")


        indices_to_flip = []
        for variant_indice, snp_name in enumerate(self.variant_names):
            own_alleles = [self.variant_dict[snp_name].minor_allele,
                              self.variant_dict[snp_name].major_allele]

            other_alleles = [other_bgen_file.variant_dict[snp_name].minor_allele,
                             other_bgen_file.variant_dict[snp_name].major_allele]


            if set(own_alleles) != set(other_alleles):
                raise ValueError(f"Alleles sets are not the same for SNP "
                                 f"{self.variant_dict[snp_name]}")

            #This relies on the assumption that the variants are ordered in the same way across the two objects
            if own_alleles != other_alleles:
                indices_to_flip.append(variant_indice)

        print(f"Flipping the alleles of {len(indices_to_flip)} variant")

        for indice in indices_to_flip:
            snp_name = other_bgen_file.variant_names[indice]
            old_major = other_bgen_file.variant_dict[snp_name].major_allele
            other_bgen_file.variant_dict[snp_name].major_allele = other_bgen_file.variant_dict[snp_name].minor_allele
            other_bgen_file.variant_dict[snp_name].minor_allele = old_major
            other_bgen_file.geno_expectations[other_bgen_file.geno_expectations[:, indice] != 3, indice] =  (
                2 - other_plink_file.genotypes[other_plink_file.genotypes[:, indice] != 3, indice]
                )

        return other_plink_file, other_plink_file.genotypes


