import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from collections import Counter, OrderedDict
import math

class PredictTemperature:
    def __init__(self, fill_mask, kmer:int, outdir:str=None):
        self.fill_mask=fill_mask
        self.kmer=kmer
        self.outdir = outdir
        self.aa, self.pred, self.pred_temp = '', [], None

    def kmer_predict(self, aa:str):
        self.aa = aa
        self.mask_seq = [(f"{self.aa} is <mask>", 1, len(aa))]
        aa_len = len(self.aa)
        if aa_len > self.kmer:
            for i in range(0, aa_len-self.kmer):
                sub_aa = self.aa[i:i+self.kmer]
                self.mask_seq.append((f"{sub_aa} is <mask>", i+1, i+self.kmer+1))
        self.pred = [self.fill_mask(_seq) for _seq, start, end in self.mask_seq]
    
    def get_temperature(self):
        scores, temp = [], []
        for sub_pred in self.pred:
            scores.append(sub_pred[0]['score'])
            try :
                t = int(sub_pred[0]['token_str'])
                temp.append(t)
            except Exception as e:
                temp.append(None)
        self.pred_temp = pd.DataFrame({'scores':scores, 'temperature':temp}, 
            dtype=np.float32)
    
    def draw_temp(self, temp:int=None, acc:str=None):
        try:
            fig, axs = plt.subplots(2, figsize=(9,4))
            fig.suptitle(f'length of protein is {len(self.aa)}, optimum Temperature is {temp}')
            x=range(len(self.pred_temp))
            axs[0].scatter(x, self.pred_temp.temperature, s=2)
            axs[0].set_ylabel('Predicted Temperature')
            axs[1].scatter(x, self.pred_temp.scores, s=2)
            axs[1].set_ylabel('Predition Scores')
            # hide x labels and tick labels for all but bottom subplot
            for ax in axs:
                ax.label_outer()
            plt.xlabel(f'Windows of k-mer sequences, kmer={self.kmer}')
            if acc and self.outdir:
                plt.savefig(os.path.join(self.outdir, f"{acc}.png"))
            else:
                plt.show()
        except Exception as e:
            print(e)
    
    def get_region(self, acc:str=None):
        seq, start, end, temp = [], [], [], []
        for (_seq, _start, _end), sub_pred in zip(self.mask_seq[1:], self.pred[1:]):
            _seq = _seq.replace(' is <mask>', '')
            try :
                t = int(sub_pred[0]['token_str'])
            except Exception as e:
                t = None
            if t:
                if len(temp) > 0:
                    if t == temp[-1]:
                        seq[-1] += _seq[-1]
                        end[-1] = _end
                    else:
                        seq.append(_seq)
                        start.append(_start)
                        end.append(_end)
                        temp.append(t)
                else:
                    seq = [_seq]
                    start = [_start]
                    end = [_end]
                    temp = [t]
        if temp:
            df = pd.DataFrame({'temperature':temp, 'start':start, 'end':end, 'seq': seq})
            if acc and self.outdir:
                df.to_csv(os.path.join(self.outdir, f"{acc}.txt"), sep="\t", index=False)
            return df
    
    def count_aa(self):
        counts = {}
        for file_name in os.listdir(self.outdir):
            if file_name.endswith('txt'):
                infile = os.path.join(self.outdir, file_name)
                df = pd.read_csv(infile, sep='\t')
                groups = df.groupby(['temperature'])
                for name, sub_df in groups:
                    temp = name[0]
                    if temp not in counts:
                        counts[temp]= Counter()
                    seq_str = ''.join(list(sub_df.seq))
                    counts[temp].update(seq_str)
        return counts
    
    def analyze_aa(self, counts, aa_pool):
        count = {}
        for aa in aa_pool:
            temp, percent = [], []
            for _temp, items in counts.items():
                temp.append(_temp)
                _count = items.get(aa, 0)
                _total = sum(items.values())
                percent.append(_count/_total)
            res = {'temperature':temp, 'percentage':percent}
            count[aa] = pd.DataFrame(res, dtype=np.float32)
        return count

