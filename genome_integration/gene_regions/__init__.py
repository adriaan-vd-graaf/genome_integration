"""
These functions are being used to prioritize genes in certain already associated regions.
"""

__author__      = "Adriaan van der Graaf"
__copyright__   = "Copyright 2017, Adriaan van der Graaf"


class StartEndRegions:
    def __init__(self, list_of_lists):
        self.gene_regions = [StartEndRegion] * len(list_of_lists)
        i = 0
        for region in list_of_lists:
            tmp = StartEndRegion(region)
            self.gene_regions[i] = tmp
            i+=1

    def in_gene_regions(self, chr, position):
        for i in self.gene_regions:
            if i.snp_in_region(chr, position):
                return True
        return False


class StartEndRegion:
    def __init__(self, list):
        self.chromosome = list[0]
        self.start = int(list[1])
        self.end = int(list[2])

    def snp_in_region(self, chr,position):
        return (self.chromosome == chr) and (self.start <= position <= self.end)
