"""
These functions are being used to help me traverse the eqtlgen things.
"""

from .. import file_utils


__author__      = "Adriaan van der Graaf"
__copyright__   = "Copyright 2017, Adriaan van der Graaf"

from . eqtl_mapping_utils import *

class EQTLGenProbeMapping:
    def __init__(self, file_name):
        self.file_name = file_name
        self.probe_to_transcript = {}
        self.mapping_list = file_utils.read_newline_separated_gz_file_into_list(file_name)[1:]
        for line in self.mapping_list:
            tmp = EQTLGenProbeTranscript(line)
            self.probe_to_transcript[tmp.probe_num] = tmp

    def probe_to_ensembl(self, probe):
        return self.probe_to_transcript[probe].ensembl_transcipt


class EQTLGenProbeTranscript:
    def __init__(self, line):
        self.split = line.split()
        self.probe_num = self.split[0]
        self.ensembl_transcipt = self.split[3]


class EnsembleMapping:
    def __init__(self, file_name):
        self.file_name = file_name
        self.transcript_to_id = {}
        self.gene_to_id = {}
        self.mapping_list = file_utils.read_newline_separated_gz_file_into_list(file_name)[1:]
        for line in self.mapping_list:
            tmp = EnsembleMappingLine(line)
            self.transcript_to_id[tmp.transcript] = tmp
            self.gene_to_id[tmp.gene] =tmp

    def ensembl_to_id(self, transcript):
        try:
            return self.transcript_to_id[transcript].ID
        except:
            return self.gene_to_id[transcript].ID

class EnsembleMappingLine:
    def __init__(self,line):
        self.split = line.split()
        self.transcript = self.split[1]
        self.gene = self.split[0]
        self.ID = self.split[2]