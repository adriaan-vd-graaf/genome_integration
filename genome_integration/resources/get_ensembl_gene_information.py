from genome_integration import gene_regions
import gzip

"""
These classes are intended to easily access and query ensembl genes information, 
But they have other uses as well, so it is possible that these will be joined with the ensembl


"""


class EnsemblGene(gene_regions.StartEndRegion):
    """
    Contains all the standard fields for gene information from ensembl.

    Attributes
    ----------
    ensg_id: str
        ensembl id

    gene_name: str
        gene name

    strand: str
        strand of the gene either `+` or `-`

    gc_percent: str
        percentage of GC bases.

    gene_type: str
        gene type

    ensembl_version:
        ensembl version of the gene.

    Methods
    -------

    None

    """

    def __init__(self,
                 ensg_id,
                 gene_name,
                 chromosome,
                 start,
                 end,
                 strand,
                 gc_percent,
                 gene_type,
                 ensembl_version):
        super().__init__([chromosome, start, end])

        self.ensg_id = ensg_id
        self.gene_name = gene_name
        self.strand = strand
        self.gc_percent = gc_percent
        self.gene_type = gene_type
        self.ensembl_version = ensembl_version

    def __repr__(self):
        return f"EnsemblGene object: {self.ensg_id}, {self.gene_name}, {self.chromosome}:{self.start}-{self.end},{self.strand}"

    def __str__(self):
        return f"{self.ensg_id}, {self.gene_name}, {self.chromosome}:{self.start}-{self.end},{self.strand}"


class EnsemblGenes:
    """
    EnsemblGene information

    This class contains many genes likely with ensembl data.

    Attributes
    ----------

    list_of_genes: list
        list of EnsemblGenes objects

    self.ensg_ids: set of strings
        Ensembl gene ids in the set.

    self.gene_names: list
        names of the genes

    self.ensg_to_full: dict
        dict of ensembl ids as keys and associated EnsemblGene information as values

    self.ensg_to_gene: dict
        dict of ensembl ids as keys and

    self.gene_to_full: dict
        gene to full


    self.genes_warned_about : set
        genes that are duplicate.

    self.allow_patch_overlapping_gene_names: bool
        allows for overlapping gene names to be available.

    Methods
    -------
    add_gene(ensembl_gene):
        add an ensembl_gene object to self.

    get_sorted_genes(self)
        return sorted genes.

    return_overlapping_regions(self, gene_region_to_check):
        return the overlapping regions of a StartEndRegion object.

    list_to_full(self, list, fail_on_bad_id=True)
        returns an ensemblGenes object with only the genes in the files.

    str_to_full(self, str)
        return the Ensembl genes object associated with a certain string.

    str_to_gene(self, str):
        return the gene name associated with a certain string.

    str_to_ensg(self, str):
        return the ensembl id associated with a certain string.





    """

    def __init__(self, allow_patch_overlapping_gene_names=False):

        self.list_of_genes = []

        self.ensg_ids = set()
        self.gene_names = set()

        self.ensg_to_full = {}

        self.gene_to_full = {}

        self.genes_warned_about = set()

        self.allow_patch_overlapping_gene_names = allow_patch_overlapping_gene_names

    def add_gene(self, ensembl_gene):
        """
        Add an EnsemblGene object to self.

        :param ensembl_gene:
        :return: None
        """
        self.list_of_genes.append(ensembl_gene)
        if ensembl_gene.ensg_id in self.ensg_ids:
            raise ValueError("ERROR: found duplicate ENSG ID, when adding {}, this should not happen.".format(ensembl_gene.ensg_id))
        self.ensg_ids.add(ensembl_gene.ensg_id)

        if ensembl_gene.gene_name in self.gene_names and ensembl_gene.gene_name not in self.genes_warned_about:
            # print("WARNING: found duplicate gene name, when adding {}, lookups on gene name may be wrong.".format(ensembl_gene.gene_name))
            self.genes_warned_about.add(ensembl_gene.gene_name)


        self.ensg_to_full[ensembl_gene.ensg_id] = ensembl_gene

        #this ensures that there will never be weird patch genes in the
        if ensembl_gene.gene_name in self.gene_names and (not self.allow_patch_overlapping_gene_names):
            try:
                len(ensembl_gene.chromosome < 3) #only chromosome names smaller than 3.
            except:
                return

        self.gene_names.add(ensembl_gene.gene_name)
        self.gene_to_full[ensembl_gene.gene_name] = ensembl_gene



    def get_sorted_genes(self):
        return sorted(self.list_of_genes)

    def return_overlapping_regions(self, gene_region_to_check):
        """
        This may be a bit slow, as it will iterate over all gene regions here.

        :param gene_region_to_check:
        :return:
        """
        sorted_genes = self.get_sorted_genes()
        to_return = EnsemblGenes()
        for gene in sorted_genes:
            if gene_region_to_check.region_overlaps(gene):
                to_return.add_gene(gene)

        return to_return


    def return_overlapping_regions_based_on_coordinates(self, chromosome, position):
        """

        This may be a bit slow, as it will iterate over all gene regions here.
        Cool thing though, this is sorted.


        :param gene_region_to_check:
        :return: EnsemblGenes object with overlapping genes.

        """
        sorted_genes = self.get_sorted_genes()
        to_return = EnsemblGenes()
        for gene in sorted_genes:
            if gene.snp_in_region(chromosome, position):
                to_return.add_gene(gene)

        return to_return


    def __str__(self):
        return "EnsemblGenes object containing {} genes".format(len(self.gene_names))

    def list_to_full(self, list, fail_on_bad_id=True):
        """
        turn a list of gene identifiers into a list of ensembl gene information

        :param list: list of IDs you want to know all the ensembl information of.
        :param fail_on_bad_id: bool,
            if bad IDs should fail. Default is True.
        :return: list of ensembl genes informaiton.
        """
        if fail_on_bad_id:
            return [self.str_to_full(x) for x in list]
        else:
            return_list = []
            for gene in list:
                try:
                    return_list.append(self.str_to_full(gene))
                except ValueError:
                    print(f"Could not find {gene}, but continueing.")

            return return_list

    def str_to_full(self, str):
        if str in self.ensg_to_full.keys():
            return self.ensg_to_full[str]
        elif str in self.gene_to_full.keys():
            return self.gene_to_full[str].gene_name
        else:
            raise ValueError(f"{str} was not convertible to a gene that I know.")

    def str_to_gene(self, str):
        return self.str_to_full(str).gene_name

    def str_to_ensg(self, str):
        return self.str_to_full(str).ensg_id

    def __iter__(self):
        self.ordered_ensembl_info = sorted(self.list_of_genes)
        self.ordered_ensg_ids = [x.ensg_id for x in self.ordered_ensembl_info]
        self.ordered_gene_names = [x.gene_name for x in self.ordered_ensembl_info]
        self.iterator_indice = 0
        return self

    def __next__(self):

        if self.iterator_indice < len(self.ordered_ensg_ids):
            self.iterator_indice += 1
            return self.ensg_to_full[self.ordered_ensg_ids[self.iterator_indice - 1]]
        else:
            raise StopIteration()



def read_gene_information():
    """
    This loads in the ENSG gene information from the package and returns it.
    very handy to have if you want to do a quick check of a certain ENSG ID, or just want gene names of everything.

    TODO: Properly handle the Ensembl and human genome versions.

    :return:
    EnsemblGene object with all the genes that are in the file '2018_05_18_ensembl_gene_information.txt.gz' in the
    resource/ensembldata folder of this package.
    """

    resource_path = '/'.join(('ensembl_data', '2018_05_18_ensembl_gene_information.txt.gz'))

    gene_file =  "{}/{}".format("/".join(__file__.split("/")[:-1]), resource_path)

    ensembl_genes = EnsemblGenes()
    with gzip.open(gene_file, "rb") as f:
        f.readline()
        for line in f:
            split = line.decode("utf8").split()
            ensembl_genes.add_gene(
                EnsemblGene(
                            split[0],
                            split[1],
                            split[2],
                            split[3],
                            split[4],
                            split[5],
                            split[6],
                            split[7],
                            split[8]
                )
            )

    return ensembl_genes