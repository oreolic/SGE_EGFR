#%%
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime

class General:
    def rt_seq(self,seq):
        rt = ''
        dic = {'A':'T','G':'C','T':'A','C':'G','N':'N','a':'t','g':'c','t':'a','c':"g"}
        for i in seq[::-1]:
            rt += dic[i]
        return rt


class DataParsing:
    def dfdic(self,df):
        dflst = self.divide_dataframe(df)
        dic = {}
        for eachdf in dflst:
            sb = eachdf.iloc[0,0]
            dic[sb] = eachdf
        return dic

    def divide_dataframe(self,df):
        bclst = self._barcode_indexing(df)
        
        t1 = datetime.now()
        n = 0
        lst = []
        while n < len(bclst)-1:
            tup0 = bclst[n]
            tup1 = bclst[n+1]
            
            idx = tup0[1]
            idx1 = tup1[1]

            #tup = (BC, idx)

            eachdf = df[idx:idx1]
            lst.append(eachdf)

            if n == len(bclst)-2:
                lastdf = df[idx1:]              
                lst.append(lastdf)

            n += 1
        t2 = datetime.now()

        print('Make DF list:', t2-t1)
        return lst

    def _barcode_indexing(self,df):
        first_column = list(df.columns)[0]
        df = df.sort_values(by=first_column)
    
        ## indexing first colunmns 
        t1 = datetime.now()
        n = 1
        
        sb = list(df.iloc[:,0])
        lst = [(sb[0],0)]
        while n < df.shape[0]:
            if sb[n] != sb[n-1]:
                lst.append((sb[n],n))
            else:
                pass
            n+=1
        t2 = datetime.now()

        print('Barcode_Indexing:', t2-t1)
        return lst
    
    def dflst(self,df):
        bclst = self._barcode_indexing(df)
        n = 0 
        t1 = datetime.now()
        n = 0
        lst = []
        while n < len(bclst)-1:
            tup0 = bclst[n]   ## tup = (BC, idx)
            tup1 = bclst[n+1]
            
            idx = tup0[1]
            idx1 = tup1[1]

            eachdf = df[idx:idx1]
            lst.append(eachdf)

            if n == len(bclst)-2:
                lastdf = df[idx1:]              
                lst.append(lastdf)

            n += 1
        t2 = datetime.now()

        print('Make DF list:', t2-t1)
        return lst

    def divide_list(self,lst, num):
        chunk_size = len(lst)//num
        final = []
        for i in range(num):
            if i == num-1:
                eachlst = lst[chunk_size*(i):]
            else:   
                eachlst = lst[chunk_size*(i):chunk_size*(i+1)]
            final.append(eachlst)
        return final



class Codon:
    def __init__(self):
        self.codon = {'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu', 'CTT': 'Leu', 'CTC': 'Leu','CTA': 'Leu',
            'CTG': 'Leu', 'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met', 'GTT': 'Val',
            'GTC': 'Val','GTA': 'Val', 'GTG': 'Val', 'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
            'CCT': 'Pro','CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro', 'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr',
            'ACG': 'Thr','GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala', 'TAT': 'Tyr', 'TAC': 'Tyr',
            'TAA': 'STOP','TAG': 'STOP', 'CAT': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln', 'AAT': 'Asn',
            'AAC': 'Asn','AAA': 'Lys', 'AAG': 'Lys', 'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
            'TGT': 'Cys','TGC': 'Cys', 'TGA': 'STOP', 'TGG': 'Trp', 'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg',
            'CGG': 'Arg','AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg', 'GGT': 'Gly', 'GGC': 'Gly',
            'GGA': 'Gly', 'GGG': 'Gly'}
        self.codon_short = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 
                        'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V','STOP':'*'}
        self.codon_abb = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'AGA': 'R', 'AGC': 'S', 'AGG': 'R',
                       'AGT': 'S', 'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I', 'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P',
                       'CCT': 'P', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E',
                       'GAT': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V',
                       'GTT': 'V', 'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TGA': '*', 'TGC': 'C', 'TGG': 'W',
                       'TGT': 'C', 'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}
        self.nucleotide = {'A':['A'],'C':['C'],'G':['G'],'R':['A','G'],'Y':['C','T'],'S':['G','C'],'W':['A','T'],'K':['G','T'],'M':['A','C'],'B':['C','G','T'],'D':['A','G','T'],'H':['A','C','T'],'V':['A','C','G'],'N':['A','T','C','G'],'T':['T']}



class Fastq:
    def read_fastq(self,path):
        dic = {}
        with open(path,'r') as fastq:
            n = 0
            while True:
                read = fastq.readline()
                if not read:
                    break
                else:
                    if n % 4 == 1:
                        read = read.rstrip()
                        if read not in dic:
                            dic[read] = 1
                        else:
                            dic[read] += 1
                    else:
                        pass
                n += 1
        return dic

    def read_fastq_example(self,path):
            dic = {}
            with open(path,'r') as fastq:
                n = 0
                while n<300000:
                    read = fastq.readline()
                    if not read:
                        break
                    else:
                        if n % 4 == 1:
                            read = read.rstrip()
                            if read not in dic:
                                dic[read] = 1
                            else:
                                dic[read] += 1
                        else:
                            pass
                    n += 1
            return dic
# %%
