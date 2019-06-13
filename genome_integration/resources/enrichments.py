import requests
import json
import sys
from genome_integration.resources.get_ensembl_gene_information import *

class Enrich:
    """
    Enrich is a class that uses the enrichr API to make enrichments of gene sets.

    Makes requests to an outside server. So don't use it if you don't want that.

    Attributes
    ----------
    list_to_enrich: list of str
        List of strings with genes to enrich.

    gene_info: EnsemblGenes
        will be filled with gene information for reference

    gene_name_list: list of str
        List of gene names after conversion.

    enrichment_json: str in json format
        json that contains the names of the enrichments.

    background: list of str
        backgrounds behind which to enrich.

    Methods
    -------
    enrich(self, enrichment_name = "None", output=True)
        Does an enrichment on  the initialized gene list

    reset(self)
        removes all the data from the class.

    convert_gene_list
        converts a gene list to names that enrichr accepts.

    do_enrichr_analysis(self, target_name)
        enriches the gene list.

    write_fdr_one_in_twenty_go_enrichments(self, file_name)
        writes the fdr < 0.05 results to a stdout (if file_name == None) or a file

    """
    def __init__(self, gene_list):
        self.list_to_enrich = list(gene_list)
        self.gene_info = None
        self.gene_name_list = None

        self.enrichment_json = None
        self.background = ['GO_Biological_Process_2017', 'Reactome_2016']

    def enrich(self, enrichment_name = "None", output=True):
        """
        Enriches the gene list.

        :param enrichment_name:
        :param output: if it should output the file.
        :return: None
        """
        if self.enrichment_json is None:
            self.convert_gene_list()
            self.do_enrichr_analysis(enrichment_name)

        if output:
            self. write_fdr_one_in_twenty_go_enrichments()

    def reset(self):

        """
        Resets to a zero gene state.

        :return: None
        """
        self.gene_info = None
        self.gene_name_list = None

        self.enrichment_json = None
        self.background = ['GO_Biological_Process_2017', 'Reactome_2016']

    def convert_gene_list(self):
        """
        internal function to return the gene list.
        :return:None
        """

        if self.gene_info is None:
            self.gene_info = read_gene_information()

        self.gene_name_list = []

        for gene in self.list_to_enrich:
            try:
                self.gene_name_list.append( self.gene_info.str_to_gene(gene))
            except ValueError as x:
                print(f"{x}, occured, continueing with the following gene" )


    def do_enrichr_analysis(self, target_name):
        """
        Does the request to enrichr for the enrichment analysis.

        :param target_name:
        :return:
        """

        ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'

        genes_str = '\n'.join(self.gene_name_list)
        description = target_name
        payload = {
            'list': (None, genes_str),
            'description': (None, description)
        }

        response = requests.post(ENRICHR_URL, files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list', target_name)

        data = json.loads(response.text)

        """
        Data is uploaded, now we request the Kegg pathways.
        """
        ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
        user_list_id = data["userListId"]

        backgrounds = self.background

        return_dict = {}
        for background in backgrounds:
            try:
                response = requests.get(
                    ENRICHR_URL + "?userListId={}&backgroundType={}".format(user_list_id, background)
                )
                if not response.ok:
                    raise Exception("Did not get correct request", target_name)

                data = json.loads(response.text)
                return_dict.update(data)
            except Exception as x:
                print("Exception raised when requesting {} as a background:".format(background))
                print(x)

        self.enrichment_json = return_dict

        return return_dict


    def write_fdr_one_in_twenty_go_enrichments(self, filename = None):

        """
        Writes the enrichr results with an FDR > 0.05 to a file or standard out.
        :param filename: if None, will write to stdout, otherwise will write to a file.
        :return: None.
        """

        if filename is not None:
            file_handle = open(filename, "w")
        else:
            file_handle = None

        print(
            "background\tTerm name\tP-value\tZ-score\tCombined_score\tOverlapping_genes\tAdjusted_p-value\n",
            file=file_handle, end=""
        )

        for key in self.enrichment_json.keys():
            all_enrichments = self.enrichment_json[key]
            for enrichment in all_enrichments:
                if enrichment[6] < 0.05:
                    string = "\t".join([str(x) for x in enrichment[1:]])
                    print(f"{key}\t{string}\n", end="", file=file_handle)

        if filename is not None:
            file_handle.close()

