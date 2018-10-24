import requests
import json
import sys
from genome_integration.resources.get_ensembl_gene_information import *

class Enrich:
    def __init__(self, gene_list):
        self.list_to_enrich = list(gene_list)
        self.gene_info = None
        self.gene_name_list = None

        self.enrichment_json = None
        self.background = ['GO_Biological_Process_2017', 'Reactome_2016']

    def enrich(self, enrichment_name = "None", output=True):

        if self.enrichment_json is None:
            self.convert_gene_list()
            self.do_enrichr_analysis(enrichment_name)

        if output:
            self.write_fdr_one_in_twenty_go_enrichments()

    def reset(self):
        self.gene_info = None
        self.gene_name_list = None

        self.enrichment_json = None
        self.background = ['GO_Biological_Process_2017', 'Reactome_2016']

    def convert_gene_list(self):

        if self.gene_info is None:
            self.gene_info = read_gene_information()

        self.gene_name_list = []

        for gene in self.list_to_enrich:
            try:
                self.gene_name_list.append( self.gene_info.str_to_gene(gene))
            except ValueError as x:
                print(f"{x}, occured, continueing with the following gene" )



    def do_enrichr_analysis(self, target_name):
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
        #bit of hack to get it to work.
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

