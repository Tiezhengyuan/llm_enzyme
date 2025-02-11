import os
import numpy as np
import pandas as pd
from Bio import SeqIO
import re
import urllib.request
import lxml.html
import json
import subprocess
import gzip


from utils import Utils

class ProcessData:
    def __init__(self, data_dir):
        self.curr_dir = data_dir
        self.meta = []

    def get_meta(self):
        return self.meta
    
    def genome_report_atcc(self, infile:str, nrows:int=None) -> list:
        '''
        get genome report from the file downloaded from NCBI
        update self.meta
        '''
        df = pd.read_csv(infile, sep='\t', na_values='', low_memory=False, nrows=nrows)
        # filter
        atcc_df = df[df.Strain.astype(str).str.contains("ATCC")]
        atcc_no = [re.findall(r'ATCC.(\d+)', i) for i in list(atcc_df.Strain)]
        atcc_df = atcc_df.assign(ATCC_No= [i[0] if i else '' for i in atcc_no])
        # convert to dict
        self.meta = atcc_df.to_dict(orient='records')
        print(f"Number of genome: {len(self.meta)}")

    def retrieve_atcc(self, atcc_dir:str):
        '''
        scan HTML downloaded from ATCC, update incubation temperature
        update self.meta
        '''
        m=n = 0
        for item in self.meta:
            strain = item.get('Strain', '')
            atcc_strain = re.findall(r'ATCC.(\d+)', strain)
            if atcc_strain and 'incubation_temperature' not in item:
                atcc_no = atcc_strain[0]
                try:
                    # download data from ATCC
                    xmlstring = self.get_atcc_data(atcc_dir, atcc_no)
                    # retrieve temperature
                    res = self.get_temperature(xmlstring)
                    item.update(res)
                    n += 1
                except Exception as e:
                    m += 1
                    # print(f"Failure: {atcc_no}, error={e}")
        print(f"parse strains={n}, failed in download={m}")


    def collect_meta(self):
        '''
        combine all meta data
        remove duplicates
        '''
        meta = []
        for name in ['genome_report_atcc.json', 'genome_report_2.json']:
            json_file = os.path.join(self.curr_dir, name)
            if os.path.isfile(json_file):
                with open(json_file, 'r') as f:
                    meta += json.load(f)
        else:
            print(f"{len(meta)}")

        # get temperature and faa
        res = {}
        for item in meta:
            # print(item)
            if 'incubation_temperature' in item:
                temperature = [str(i) for i in item['incubation_temperature']]
                for name in os.listdir(item.get('local_path')):
                    if name.endswith('.faa.gz'):
                        key = item['#Organism/Name']
                        faa_file = os.path.join(item.get('local_path'), name)
                        res[key] = (key, temperature, faa_file)
        return res

    def prepare_aa(self, aa_meta):
        '''
        '''
        m, n = 0, 0
        # create directories
        outdir = os.path.join(self.curr_dir, 'temperature')
        for specie, temperature, faa_file in aa_meta.values():
            temperature_str = '-'.join(temperature)
            m += 1
            outdir2 = os.path.join(outdir, str(temperature))
            if not os.path.isdir(outdir2):
                os.mkdir(outdir2)
            # create data
            outfile = os.path.join(outdir2, specie)
            with open(outfile, 'w') as f:
                x = 0
                for seq_id, seq_desc, seq in self.read_faa(faa_file):
                    _item = [seq, temperature_str, specie, seq_id, seq_desc, '\n']
                    f.write(' '.join(_item))
                    x += 1
                if x==0:
                    print(f"Warning: not enzyme is detected in {faa_file}")
                #     break
                n += x
        print(f"species={m}, proteins={n}")

    def prepare_aa_mask(self, aa_meta):
        m, n = 0, 0
        # create directories
        outdir = os.path.join(self.curr_dir, 'temperature_mask')
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        outfile = os.path.join(outdir, f"temperature_mask.txt")
        with open(outfile, 'w') as f:
            for specie, temperature, faa_file in aa_meta.values():
                temperature = np.max([int(t) for t in temperature])
                m += 1
                for seq_id, seq_desc, seq in self.read_faa(faa_file):
                    _item = [seq, ' is ', str(temperature), '\n']
                    f.write(' '.join(_item))
                    n += 1
        print(f"species={m}, proteins={n}")
                

    def read_faa(self, faa_file):
        with gzip.open(faa_file, 'rt') as f:
            for rec in SeqIO.parse(f, 'fasta'):
                if 'ase' in rec.description:
                    yield rec.id, rec.description, str(rec.seq)


    
    def get_atcc_data(self, indir:str, atcc_no):
        '''
        read HTMl. Download it if that doesn't exist.
        '''
        html_file = os.path.join(indir, f'{atcc_no}.html')
        if os.path.isfile(html_file):
            with open(html_file, 'r') as f:
                return f.read()

        # download
        user_agent = 'Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.9.0.7) Gecko/2009021910 Firefox/3.0.7'
        url = f"https://www.atcc.org/products/{atcc_no}"
        headers={'User-Agent':user_agent,} 

        #The assembled request
        request=urllib.request.Request(url,None,headers) 
        response = urllib.request.urlopen(request)
        data = response.read()
        xmlstring = data.decode('utf-8')
        # 
        with open(html_file, 'w') as f:
            f.write(xmlstring)
        return xmlstring
    
    def get_temperature(self, xmlstring):
        res = {}
        html = lxml.html.fromstring(xmlstring)
        for el in html.iter():
            if el.text and '°C' in el.text and el.tag=='div':
                _temp = [int(i) for i in re.findall(r'(\d+)', el.text)]
                if 'to' in el.text:
                    res['storage_temperature'] = _temp
                else:
                    res['incubation_temperature'] = _temp
        return res


