import os
import json
from Bio import SeqIO
import gzip

class Utils:

    @staticmethod
    def from_json(infile:str):
        try:
            with open(infile, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(e)
        return {}

    @staticmethod
    def to_json(data:dict, outdir:str, file_name:str):
        json_file = os.path.join(outdir, file_name)
        with open(json_file, 'w') as f:
            json.dump(data, f, indent=True)
        return json_file
    
    @staticmethod
    def iter_uniprot(gz_file:str):
        if not os.path.isfile(gz_file):
            return None

        with gzip.open(gz_file, 'rt') as f:
            for rec in SeqIO.parse(f, 'swiss'):
                yield rec