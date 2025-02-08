
import os
import pandas as pd
from Bio import SeqIO
import re
import urllib.request
import lxml.html
import json
import subprocess

def download_genome(ftp_path):
    local_dir = '/home/yuan/data'
    try:
        subprocess.run(['wget', '-r', '-P', local_dir, ftp_path])
        _url = ftp_path.replace('ftp://', '')
        return os.path.join(local_dir, _url)
    except Exception as e:
        print(e)


if __name__ == "__main__":
    outfile = 'data/genome_report_2.json'
    infile = outfile if os.path.isfile(outfile) else 'data/genome_report.json'
    with open(infile, 'r') as f:
        data = json.load(f)
    
    for item in data:
        if 'local_path' not in item:
            print(item['FTP Path'])
            local_path = download_genome(item['FTP Path'])
            if os.path.isdir(local_path):
                item['local_path'] = local_path
                print(item['local_path'])
        else:
            print(f"Skip {item['FTP Path']}")

    with open(outfile, 'w') as f:
        json.dump(data, f, indent=True)
