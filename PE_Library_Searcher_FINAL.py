#%%
import pandas as pd
from concurrent.futures.process import ProcessPoolExecutor
from Codes import general as gen
from Codes import umi_collapse as ucs
from datetime import datetime
import os

#%%

def edit_dist(str1, str2):
    str1 = str1.upper()
    str2 = str2.upper()
    dp = [[0] * (len(str2)+1) for _ in range(len(str1) + 1)]
    for i in range(1, len(str1)+1):
        dp[i][0] = i
    for j in range(1, len(str2)+1):
        dp[0][j] = j

    for i in range(1, len(str1)+1):
        for j in range(1, len(str2)+1):
            if str1[i-1] == str2[j-1]:
                dp[i][j] = dp[i-1][j-1]

            else:
                dp[i][j] = min(dp[i-1][j-1], dp[i-1][j], dp[i][j-1]) + 1

    return dp[-1][-1]


def _ed(df):
    rt = df['RTPBS'].upper()
    ref = df['REF'].upper()

    edit = edit_dist(rt,ref)
    return edit

def ed(eachdf):
    t1 = datetime.now()
    eachdf['ED'] = eachdf.apply(_ed,axis=1)
    t2 = datetime.now()
    print('Calculate ED  ',t2-t1)
    return eachdf


class ReferencePosition:
    def __init__(self):
        self.ref_scaffold = -1
        self.ref_barcode = -1
        self.ref_rp = -1
        self.ref_UMI = -1
        self.barcode_length = -1
        self.umi_length = -1
        self.scaffold_length = -1
        self.rp_length = -1
        self.rtpbs_index = -1

        return 

def defining_reference_position(rs):
    inp = pd.read_csv('Input/library_analyzer_input.txt',sep='\t')
    try:
        ref = inp.loc[0,'Reference']
        scaf = inp.loc[0,'Scaffold']
        barcode = inp.loc[0,'Barcode']
        umi = inp.loc[0,'UMI']
        rp = inp.loc[0,'RP']
        rtpbs = inp.loc[0,'RTPBS']
    except:
        print('Check Input Files, Column names sholud be Reference/Scaffold/Barcode/UMI/RP, and Input sequences')
        return 'FALSE'
    
    scaf_index = ref.find(scaf)
    barcode_index = ref.find(barcode)
    umi_index = ref.find(umi)
    rp_index = ref.find(rp)


    rs.ref_scaffold = scaf_index
    rs.ref_barcode = barcode_index
    rs.ref_rp =rp_index
    rs.umi_index = umi_index
    rs.barcode_length = len(barcode)
    rs.umi_length = len(umi)
    rs.scaffold_length = len(scaf)
    rs.rp_length = len(rp)
    rs.scaffold_seq = scaf
    rs.rtpbs_index = ref.find(rtpbs)
    
    return rs
        
def _generate_ed1_bc(bclst):
    ed1_bc_dic = {}
    for idx in bclst:
        ed1_bc_dic[idx] = idx
        ed1lst = _generate_ed1_seq(idx)
        for ed1bc in ed1lst:
            ed1_bc_dic[ed1bc] = idx

    return ed1_bc_dic

def _generate_ed1_seq(seq):
    base = ['A', 'T', 'G', 'C']
    ed1_list = [seq]
    for i in range(len(seq)):
        if i == 0:
            for nt in base:
                ed1 = nt + seq[1:]
                if ed1 == seq:
                    pass
                else:
                    ed1_list.append(ed1)
        elif i == len(seq) - 1:
            for nt in base:
                ed1 = seq[:-1] + nt
                if ed1 == seq:
                    pass
                else:
                    ed1_list.append(ed1)
        else:
            for nt in base:
                ed1 = seq[:i] + nt + seq[i + 1:]
                if ed1 == seq:
                    pass
                else:
                    ed1_list.append(ed1)

    return ed1_list

def _find_barcode(read,ed1_bc,barcode_length):
    ## ed1_bc = dictinary {}
    n= 0
    sb = 'FALSE'
    barcode_index = -1
    while n < len(read)-barcode_length+1:
        motif = read[n:n+barcode_length]
        if motif in ed1_bc:
            sb = ed1_bc[motif]
            barcode_index = n
            break
        else: pass
        n += 1

    return sb,barcode_index

def find_barcode_in_fastqlst(eachfastqlst,ed1_bc,barcode_length):
    t1 = datetime.now()
    lst = []
    dic = {}
    for tup in eachfastqlst:
        read = tup[0]
        count = tup[1]

        sb, sb_index = _find_barcode(read,ed1_bc,barcode_length)
        if sb == 'FALSE':
            pass
        else:
            dic[read] = [sb,sb_index,count]
            lst.append([read,sb,sb_index,count])

    eachdf = pd.DataFrame(lst,columns = ['READ','Barcode','Barcode_Index','ReadCount'])
    t2 = datetime.now()
    print('Find Barcode: ',t2-t1)
    return eachdf,dic

def find_barcode(fastqlst,bc):
    ed1_bc = _generate_ed1_bc(list(bc.index))
    barcode_length = len(list(ed1_bc.keys())[0])
    num = 40
    with ProcessPoolExecutor(max_workers=num) as executor:
        futs = []
        chunk_size = len(fastqlst)//num
        for i in range(num):
            if i == num-1:
                eachfastqlst = fastqlst[chunk_size*(i):]
            else:   
                eachfastqlst = fastqlst[chunk_size*(i):chunk_size*(i+1)]

            fut = executor.submit(find_barcode_in_fastqlst,eachfastqlst,ed1_bc,barcode_length)
            futs.append(fut)

        merged_dic = {}
        for fut in futs:
            eachdf, eachdic = fut.result()
            for read in eachdic:
                merged_dic[read] = eachdic[read]
               
        return merged_dic

def _find_rp(read,ed1_rp,rp_length):
    ## ed1_bc = dictinary {}
    n= 0
    rp = 'FALSE'
    rp_index = -1
    while n < len(read)-rp_length+1:
        motif = read[n:n+rp_length]
        if motif in ed1_rp:
            rp = ed1_rp[motif]
            rp_index = n
            break
        else: pass
        n += 1

    return rp,rp_index

def find_rp_in_fastqlst(eachfastqlst,ed1_rp,rp_length):
    t1 = datetime.now()
    dic = {}
    for tup in eachfastqlst:
        read = tup[0]
        eachresult = tup[1] ## eachresult = [sb, sbindex,readcount,scaffold_index]

        rp,rp_index = _find_rp(read,ed1_rp,rp_length)
        if rp == 'FALSE':
            pass
        else:
            eachresult.extend([rp,rp_index])
            dic[read] = eachresult


    t2 = datetime.now()
    print('Find RP: ',t2-t1)
    return dic

def find_rp(fastqlst,bc):
    ed1_rp = _generate_ed1_bc(list(set(bc['RP'])))
    rp_length = 20

    num = 40
    with ProcessPoolExecutor(max_workers=num) as executor:
        futs = []
        chunk_size = len(fastqlst)//num
        for i in range(num):
            if i == num-1:
                eachfastqlst = fastqlst[chunk_size*(i):]
            else:   
                eachfastqlst = fastqlst[chunk_size*(i):chunk_size*(i+1)]

            fut = executor.submit(find_rp_in_fastqlst,eachfastqlst,ed1_rp,rp_length)
            futs.append(fut)

        merged_dic = {}
        for fut in futs:
            eachdic = fut.result()
            for read in eachdic:
                merged_dic[read] = eachdic[read]
        
        return merged_dic
    
def _define_umi(df):
    rpseq = df['RP']
    rp_length = len(rpseq)

    rp_index = df['RP_Index']
    read = df['Read']
    return read[rp_index+rp_length:rp_index+rp_length+8]

def define_rtpbs_umi(dic):
    ## dic = {read:[barcode, barcode index, read count, scaffold_index,rp,rp_index]}
    lst = []
    
    for read in dic:
        eachresult = dic[read]
        eachresult.append(read)
        lst.append(eachresult)
    
    result = pd.DataFrame(lst,columns = ['Barcode','Barcode_Index','ReadCount','RP','RP_Index','Read'])

    result['UMI'] = result.apply(_define_umi,axis=1)

    return result




#%%
def main(fastq,output,bc):
    ## Defining NGS Read Component Reference Positio
    #bc = bc.set_index('SortingBarcode')
    dic = gen.Fastq().read_fastq(fastq)
    readlst = list(dic.items())

    ### Searching Barcode
    dic = find_barcode(readlst,bc)
    readlst = list(dic.items())
    ## dic = {read:[barcode, barcode index, read count]} 

    ## Searching RP
    dic = find_rp(readlst,bc)
    ## dic = {read:[barcode, barcode index, read count, scaffold_index,rp,rp_index]}

    ### Define RTPBS and UMI
    result = define_rtpbs_umi(dic)


    result = result.sort_values(by='Barcode')
    result = result[['Barcode','UMI','ReadCount']]
    result = result.sort_values(by='Barcode')
    result.index = [ i for i in range(result.shape[0])]
    dfdic = gen.DataParsing().dfdic(result)
    
    result = ucs.umi_collapsing_for_library(dfdic)

    result.to_pickle('Output/{}_Barcode_UMIcollpased.pkl'.format(output))

    return

if __name__ == "__main__":
    bc = pd.read_csv('Barcode/EGFR_Library_Barcode_Final.txt',sep='\t',index_col=1)
    fastq = 'FASTQ/example.fastq'
    output = 'example'
    files = os.listdir('FASTQ/T790M')

    main(fastq,output,bc)

# %%
