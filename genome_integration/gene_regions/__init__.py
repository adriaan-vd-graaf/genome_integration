"""
These functions are being used to prioritize genes in certain already associated regions.

TODO: add a gene region with splicing in it.
"""

__author__      = "Adriaan van der Graaf"
__copyright__   = "Copyright 2017, Adriaan van der Graaf"


class StartEndRegions:
    def __init__(self, list_of_regions):
        self.gene_regions = [StartEndRegion] * len(list_of_regions)
        i = 0
        for region in list_of_regions:
            tmp = StartEndRegion(region)
            self.gene_regions[i] = tmp
            i+=1

    def in_gene_regions(self, chr, position):
        for i in self.gene_regions:
            if i.snp_in_region(chr, position):
                return True
        return False

    def make_non_overlapping_regions(self):

        already_combined = set()

        new_list = []

        for i in range(len(self.gene_regions)):
            if i in already_combined:
                continue

            region = self.gene_regions[i]
            overlapping = [j for j in range(len(self.gene_regions)) if self.gene_regions[j].region_overlaps(region)]

            chromosome = region.chromosome
            start = min([self.gene_regions[j].start for j in overlapping])
            end = max([self.gene_regions[j].end for j in overlapping])
            # add it to a list.
            new_list.append([chromosome, start, end])
            # finally remove the regions that do not have overlap.
            [already_combined.add(j) for j in overlapping]

        return StartEndRegions(new_list)


class StartEndRegion:
    def __init__(self, list):
        self.chromosome = str(list[0])
        self.start = int(list[1])
        self.end = int(list[2])

        # Runtime checks.
        if self.start > self.end:
            raise RuntimeError("Region cannot have a start position smaller than an end position")

        if self.start < 0 or self.end < 0:
            raise RuntimeError("Region cannot have negative positions.")

    def snp_in_region(self, chr, position):
        return (self.chromosome == str(chr)) and (self.start <= int(position) <= self.end)

    def snp_object_in_region(self, snp_object):
        """
        :param snp_object:
        :return: True or false if snp in region
        """
        return self.snp_in_region(snp_object.chromosome, snp_object.position)


    def region_overlaps(self, other_region):
        if self.chromosome == other_region.chromosome:
            #this may contain an error, and could be done more efficiently.
            if self.start <= other_region.start <= self.end \
                    or self.start <= other_region.end <= self.end\
                    or other_region.start <= self.start <= other_region.end \
                    or other_region.start <= self.end <= other_region.end:
                return True

        return False

    def __str__(self):
        return '{}:{}-{}'.format(self.chromosome, self.start, self.end)

    def __lt__(self, other):

        if not other.__class__ is self.__class__:
            return NotImplemented

        if not self.chromosome == other.chromosome:
            try:
                return int(self.chromosome) < int(other.chromosome)
            except:
                return self.chromosome < other.chromosome

        return self.start < other.start
