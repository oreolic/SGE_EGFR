import pandas as pd
import itertools
import Endo_REF
import os
from itertools import combinations
from itertools import product

def read_fastqntrim(file, ref):
    
    FP = ref[:10]
    FP_lst = [FP]
    nucleo_dic = {"A": 0,"T": 0,"G": 0,"C": 0}

    for i in range(len(FP)):  

        for key in nucleo_dic: 
            if FP[i] == key:
                pass
            else:
                nucleo_var = FP[0:i] + key + FP[i+1:]
                FP_lst.append(nucleo_var)

    input_dic = {}
    with open("{}".format(file), 'r') as f:
        for i, line in enumerate(f):
            if i%4 !=1: pass
            else:
                for FP in FP_lst:
                    if line.find(FP) != -1:
                        x = line[(line.find(FP)+len(FP)):(line.find(FP)+len(FP)+len(ref[10:]))]
                        if x in input_dic.keys(): 
                            input_dic[x] += 1
                        else:
                            input_dic[x] = 1
                    else: pass
    sub_input ={} 
    for key in input_dic:
        if len(key) == len(ref[10:]):
            sub_input[key] = input_dic[key]
        else: pass                
    

    return sub_input


def mismatch (input, n):

    input_len = len(input)
    a = list(input)
    dic = {}
    loc = list(combinations(range(input_len),n))
    nucleo_dic = {"A": ["T","G","C"], "T":["A","G","C"], "G":["A","T","C"], "C":["A","T","G"]}


    for i in loc:
        b = a.copy()
        for k in range(len(i)):
            b[i[k]] = nucleo_dic[b[i[k]]]
        lst = list(itertools.product(*b))
        for i in lst:
            dic [''.join(i)] = input
            
    return dic


def conv(R1, EE, ref):
    

    dic = {}

    for i in EE.index:
        seq = EE.loc[i,"sequence"]
        if seq not in dic:
            dic[seq] = 0
        
        else: pass

    for k in R1:
        if k in dic:
            dic[k] = R1[k]
        
        else: pass

    df = pd.DataFrame(dic.items(), columns=["sequence","ReadCounts"])



    ref_dic = {}
    aa_dic = {}
    for i in EE.index:
        seq = EE.loc[i,"sequence"]
        ID = EE.loc[i,"ID"]
        AA = EE.loc[i,"AA_change"]
        ref_dic[seq] = ID
        aa_dic[seq] = AA
    
    aa_df = pd.DataFrame(aa_dic.items(), columns=["sequence", "AA_change"])
    ref_df = pd.DataFrame(ref_dic.items(), columns=["sequence", "ID"])
    df =df.set_index("sequence")
    ref_df = ref_df.set_index("sequence")
    aa_df = aa_df.set_index("sequence")
    c = pd.concat([ref_df, aa_df, df], axis = 1)
    c = c.reset_index()

    ddic = {}
    a_dic = {}
    
    for k in c.index:
        ID = c.loc[k,"ID"]
        reads = c.loc[k,"ReadCounts"]
        aa_id = c.loc[k,"AA_change"]
        
        a_dic[ID] = aa_id
        
        if ID not in ddic:
            ddic[ID] = reads
        
        else:
            ddic[ID] += reads
    R1_WT = R1[ref[10:]]

    ddic["WT"] = R1_WT
    a_dic["WT"] = "WT"
    
    a_fin = pd.DataFrame(a_dic.items(), columns=["ID", "AA_change"]) 
    b_fin = pd.DataFrame(ddic.items(), columns=["ID", "ReadCounts"]) 

    a_fin = a_fin.set_index("ID")
    b_fin = b_fin.set_index("ID")    

    final = pd.concat([a_fin,b_fin], axis = 1)
    final = final.reset_index()
    
    for i in final.index:
        ID = final.loc[i,"ID"]
        if "_" in ID:
            site = ID.split('_')[1]
            final.loc[i,"site"] = site
        
        else: 
            final.loc[i,"site"] = "WT"


    return final


ref = Endo_REF.exon(20)
Exon = "E20"
EE = pd.read_csv("Barcode/PC9_E20_aa.csv")


R1 = read_fastqntrim("#####.fastq", ref)
R1_df = conv(R1, EE, ref)


