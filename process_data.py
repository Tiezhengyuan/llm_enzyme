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

class ProcessData:
    def __init__(self,):
        self.meta = []
        self.curr_dir = os.path.dirname(__file__)

    def get_meta(self):
        meta = []
        for name in ['genome_report_2.json', 'genome_report.json']:
            json_file = os.path.join(self.curr_dir, 'data', name)
            if os.path.isfile(json_file):
                with open(json_file, 'r') as f:
                    meta = json.load(f)
                    break
        else:
            print(f"{len(meta)}")

        # get temperature and faa
        res = []
        for item in meta:
            # print(item)
            if 'incubation_temperature' in item:
                temperature = [str(i) for i in item['incubation_temperature']]
                for name in os.listdir(item.get('local_path')):
                    if name.endswith('.faa.gz'):
                        faa_file = os.path.join(item.get('local_path'), name)
                        res.append((item['#Organism/Name'], temperature, faa_file))
        return res

    def prepare_aa(self, aa_meta):
        m, n = 0, 0
        # create directories
        outdir = os.path.join(self.curr_dir, 'data', 'temperature')
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        for specie, temperature, faa_file in aa_meta:
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
        outdir = os.path.join(self.curr_dir, 'data', 'temperature_mask')
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        outfile = os.path.join(outdir, f"temperature_mask.txt")
        with open(outfile, 'w') as f:
            for specie, temperature, faa_file in aa_meta:
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


    def retrieve_atcc(self, num_rec:int=None):
        n = 0
        for item in self.meta:
            n += 1
            strain = item.get('Strain', '')
            atcc_strain = re.findall(r'ATCC.(\d+)', strain)
            if atcc_strain and 'incubation_temperature' not in item:
                atcc_no = atcc_strain[0]
                try:
                    # download data from ATCC
                    xmlstring = self.get_atcc_data(atcc_no)
                    # retrieve temperature
                    res = self.get_temperature(xmlstring)
                    item.update(res)
                except Exception as e:
                    print(f"Failure: {atcc_no}, error={e}")
            if num_rec and n >= num_rec:
                break

        json_file = 'data/genome_report.json'
        with open(json_file, 'w') as f:
            json.dump(self.meta, f, indent=True)
        return self.meta

    def genome_report(self, infile:str, nrows:int=None) -> list:
        '''
        get genome report from the file downloaded from NCBI
        '''
        df = pd.read_csv(infile, sep='\t', na_values='', low_memory=False, nrows=nrows)
        # filter
        atcc_df = df[df.Strain.astype(str).str.contains("ATCC")]
        atcc_no = [re.findall(r'ATCC.(\d+)', i) for i in list(atcc_df.Strain)]
        atcc_df = atcc_df.assign(ATCC_No= [i[0] if i else '' for i in atcc_no])
        # convert to dict
        self.meta = atcc_df.to_dict(orient='records')
        print(f"Number of genome: {len(self.meta)}")
        return self.meta
    
    def get_atcc_data(self, atcc_no):
        html_file = f'data/atcc/{atcc_no}.html'
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


