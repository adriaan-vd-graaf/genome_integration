"""
These functions are being used to prioritize genes in certain already associated regions.
"""

__author__      = "Adriaan van der Graaf"
from .. import variants

class StartEndRegion:
    """
    Class that implements a genomic region, which is really a stretch of base pairs.

    Attributes
    ----------

    chromosome: str
        string representing the chromosome

    start: int
        start base pair position of the genomic region.

    end: int
        end position of the genomic region.

    Methods
    -------

    position_in_region(self, chr, position):
        returns a bool on if the single position is within the region.

    snp_in_region(chr, position)
        returns a bool on if the single position is within the region.
        synonym method of position in region.

    snp_object_in_region(snp_object)
        returns a bool on if the object of class SNP is within the region.

    region_overlaps(other_region)
        returns a bool if another of region of the same class overlaps.



    """
    def __init__(self, *args, **kwargs):

        if type(args[0]) == list and len(args) == 1:
            self.chromosome = str(args[0][0])
            self.start = int(args[0][1])
            self.end = int(args[0][2])
        elif len(args) == 3 and type(args[0]) != list:
            self.chromosome = str(args[0])
            self.start = int(args[1])
            self.end = int(args[2])

        elif type(args[0]) == str:
            self.chromosome, _ , rest = args[0].partition(":")
            self.start, self.end = [int(x) for x in rest.split('-')]
        else:
            raise ValueError("Constructor only accepts a list [<chr>, <start>, <end>], three arguments (<chr>, <start>, <end>) "
                             "or a string formatted as 'chr:start-end' ")

        # Runtime checks.
        if self.start > self.end:
            raise RuntimeError("Region cannot have a start position smaller than an end position")

        if self.start < 0 or self.end < 0:
            raise RuntimeError("Region cannot have negative positions.")

    def position_in_region(self, chr, position):
        """
        return if the position is in the region.

        :param chr: chromosome str or castable to str
        :param position: position int or castable to int
        :return: bool
        """
        return (self.chromosome == str(chr)) and (self.start <= int(position) <= self.end)

    def snp_in_region(self, chr, position):
        """
        return if the position is in the region.
        synonym method of position_in_region method.

        :param chr: chromosome str or castable to str
        :param position: position int or castable to int
        :return: bool
        """
        return self.position_in_region(chr, position)


    def snp_object_in_region(self, snp_object):
        """

        :param snp_object: object of type SNP (this package)
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

    def __contains__(self, item):
        if isinstance(item, variants.SNP):
            return self.snp_object_in_region(item)
        elif isinstance(item, StartEndRegion):
            return self.region_overlaps(item)
        else:
            raise ValueError("Only classes (or inheritance allowed:) SNP.variant or gene_regions.StartEndRegion")





class StartEndRegions:
    """
    This class contains multiple start end regions

    Attributes
    ----------

    gene_regions: list of StartEndRegion objects

    Methods
    -------

    in_gene_regions(self, chr, position)
        identifies if a position is in any of the regions.

    make_non_overlapping_regions(self)
        combines the regions into contiguous non-overlapping regions.

    """
    def __init__(self, list_of_regions):
        self.gene_regions = list([StartEndRegion] * len(list_of_regions))
        i = 0
        for region in list_of_regions:
            tmp = StartEndRegion(region)
            self.gene_regions[i] = tmp
            i+=1

    def in_gene_regions(self, chr, position):
        """
        Identify if a snp is in any of the gene regions.

        :param chr: chromosome castable to str
        :param position: position castable to int
        :return: boolean
        """
        for i in self.gene_regions:
            if i.snp_in_region(chr, position):
                return True
        return False

    def make_non_overlapping_regions(self):
        """
        Combines all overlapping regions, and turns them into one big region.

        :return: new instance of StartEndRegions containing contiguous, non overlapping regions
        """
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

    def __next__(self):
        self.i +=1
        if self.i > len(self.gene_regions):
            raise StopIteration
        else:
            return self.gene_regions[self.i-1]

    def __iter__(self):
        self.i = 0
        return self
