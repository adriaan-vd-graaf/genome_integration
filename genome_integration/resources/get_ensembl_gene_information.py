from genome_integration import gene_regions
import gzip


class EnsemblGene(gene_regions.StartEndRegion):
    """
    This data class (can be made into one in python 3.7, could be cool)
    contains all the standard fields for gene information from ensembl.
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

    def __str__(self):
        return(f"{self.ensg_id}, {self.gene_name}, {self.chromosome}:{self.start}-{self.end},{self.strand}")


class EnsemblGenes:
    def __init__(self, allow_patch_overlapping_gene_names=False):

        self.list_of_genes = []

        self.ensg_ids = set()
        self.gene_names = set()

        self.ensg_to_full = {}
        self.ensg_to_gene = {}

        self.gene_to_full = {}
        self.gene_to_ensg = {}

        self.genes_warned_about = set()

        self.allow_patch_overlapping_gene_names = allow_patch_overlapping_gene_names

    def add_gene(self, ensembl_gene):

        self.list_of_genes.append(ensembl_gene)
        if ensembl_gene.ensg_id in self.ensg_ids:
            raise ValueError("ERROR: found duplicate ENSG ID, when adding {}, this should not happen.".format(ensembl_gene.ensg_id))
        self.ensg_ids.add(ensembl_gene.ensg_id)

        if ensembl_gene.gene_name in self.gene_names and ensembl_gene.gene_name not in self.genes_warned_about:
            # print("WARNING: found duplicate gene name, when adding {}, lookups on gene name may be wrong.".format(ensembl_gene.gene_name))
            self.genes_warned_about.add(ensembl_gene.gene_name)


        self.ensg_to_full[ensembl_gene.ensg_id] = ensembl_gene
        self.ensg_to_gene[ensembl_gene.ensg_id] = ensembl_gene.gene_name

        #this ensures that there will never be weird patch genes in the
        if ensembl_gene.gene_name in self.gene_names and (not self.allow_patch_overlapping_gene_names):
            try:
                len(ensembl_gene.chromosome < 3) #only chromosome names smaller than 3.
            except:
                return

        self.gene_names.add(ensembl_gene.gene_name)
        self.gene_to_full[ensembl_gene.gene_name] = ensembl_gene
        self.gene_to_ensg[ensembl_gene.gene_name] = ensembl_gene.ensg_id


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

    def __str__(self):
        return "EnsemblGenes object containing {} genes".format(len(self.gene_names))

    def list_to_full(self, list, fail_on_bad_id=False):

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
        elif str in self.gene_to_ensg.keys():
            ensg = self.gene_to_ensg[str]
            return self.ensg_to_full[ensg]
        else:
            raise ValueError(f"{str} was not convertible to a gene that I know.")

    def str_to_gene(self, str):
        return self.str_to_full(str).gene_name

    def str_to_ensg(self, str):
        return self.str_to_full(str).ensg_id


    def __iter__(self):
        self.ordered_ensg_ids = list(self.ensg_ids)
        self.ordered_gene_names = list(self.gene_names)
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
    very handy to have if you want to do a quick check of a certain ENSG number.

    :return:
    EnsemblGene object with all the genes that are in the file '2018_05_18_ensembl_gene_information.txt.gz' in the
    resource/ensembldata folder.
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