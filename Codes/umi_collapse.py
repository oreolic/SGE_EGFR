#%%
from datetime import datetime
import pandas as pd
import pandas as pd
from concurrent.futures.process import ProcessPoolExecutor

       
def _generate_ed1_seq (seq):
    nts = ['A', 'T', 'G', 'C']
    if seq.find('N') != -1:
        return []
    else:
        ed1_lst = []
        for i in range(len(seq)):
            if i == 0:
                for nt in nts:
                    ed1 = nt + seq[1:]
                    if ed1 == seq:
                        pass
                    else:
                        ed1_lst.append(ed1)
            elif i == len(seq) - 1:
                for nt in nts:
                    ed1 = seq[:-1] + nt
                    if ed1 == seq:
                        pass
                    else:
                        ed1_lst.append(ed1)
            else:
                for nt in nts:
                    ed1 = seq[:i] + nt + seq[i + 1:]
                    if ed1 == seq:
                        pass
                    else:
                        ed1_lst.append(ed1)

        return ed1_lst
def _collapsing(node_branch_dic_final):

    node_branch_dic = node_branch_dic_final.copy()

    merged_dic = {}
    collaped_umi_dic = {}
    for tup in node_branch_dic:
        if tup in collaped_umi_dic: 
            pass
        else:
            brn_umi_lst = node_branch_dic[tup]
            for brn_tup in brn_umi_lst:
                if brn_tup in node_branch_dic:
                    if brn_tup in collaped_umi_dic: 
                        pass
                    else:
                        brn_umi_lst.extend(node_branch_dic[brn_tup])
                        collaped_umi_dic[brn_tup] = 1
                else: pass
            brn_umi_lst = list(set(brn_umi_lst))
            merged_dic[tup] = brn_umi_lst
    return merged_dic, collaped_umi_dic

def _generate_node_branch_dic(df):
    umi_count_dic = {}
    for idx in df.index:
        umi = df.loc[idx,'UMI']
        if umi not in umi_count_dic:
            umi_count_dic[umi] = df.loc[idx,'ReadCount']
        else:
            umi_count_dic[umi] += df.loc[idx,'ReadCount']
    umi_count_tup = list(umi_count_dic.items())
    umi_count_tup = sorted(umi_count_tup, key=lambda x:-x[1])

    node_branch_dic = {}
    for idx in range(len(umi_count_tup)):
        tup = umi_count_tup[idx]

        umi = tup[0]
        count = tup[1]

        ed1_umi = _generate_ed1_seq(umi)
        branch_umi = umi_count_tup[idx+1:]
        branch_umi = dict(branch_umi)
        branch_umi_lst = [(i, branch_umi[i]) for i in branch_umi if i in ed1_umi]
        branch_umi_lst = [br_tup for br_tup in branch_umi_lst if count > 3*br_tup[1]-1]
        node_branch_dic[tup] = list(set(branch_umi_lst))

    used_umi = {}
    node_branch_dic_final = {}
    for tup in node_branch_dic:
        used_umi[tup] = 1
        eachlst = node_branch_dic[tup]
        eachlst = [i for i in eachlst if i not in used_umi]
        for i in eachlst:
            used_umi[i] = 1
        
        node_branch_dic_final[tup] = eachlst
        
    while True:
        node_branch_dic_final,collapsed = _collapsing(node_branch_dic_final)
        if len(collapsed) == 0:
            final_node_branch_dic = node_branch_dic_final.copy()
            return final_node_branch_dic
            
def _umi_collapsing(df):
    node_branch_dic = _generate_node_branch_dic(df)
    collapsing_dic = {}
    for tup in node_branch_dic:
        collapsing_dic[tup[0]] = tup[0]
        eachlst = node_branch_dic[tup]
        for br_umi_tup in eachlst:
            br_umi = br_umi_tup[0]
            collapsing_dic[br_umi] = tup[0]

    df['UMI'] = [collapsing_dic[i] for i in df['UMI']]
    return df


def umi_collapsing_for_single_barcode(eachdf):
    eachdf = eachdf.sort_values(by='UMI')
    eachdf = eachdf.copy()

    eachdf.index = [i for i in range(eachdf.shape[0])]
    eachdf = _umi_collapsing(eachdf)
    umi_count_dic = {}
    for idx in eachdf.index:
        umi = eachdf.loc[idx,'UMI']
        count = eachdf.loc[idx,'ReadCount']
        if umi not in umi_count_dic:
            umi_count_dic[umi] = count
        else:
            umi_count_dic[umi] += count
    
    eachdf = eachdf.drop_duplicates('UMI')
    eachdf['ReadCount'] = [umi_count_dic[i] for i in eachdf['UMI']]

    return eachdf

def umi_collapsing_for_library(dfdic):
    t1 = datetime.now()
    

    with ProcessPoolExecutor(max_workers=40) as executor:
        futs = []
        for sb in dfdic:
            eachdf = dfdic[sb]
            fut = executor.submit(umi_collapsing_for_single_barcode,eachdf)
            futs.append(fut)

        merged = []
        for fut in futs:
            merged.append(fut.result())
        
        final = pd.concat(merged)
        t2 = datetime.now()
        print('UMI COLLAPSING: ', t2-t1) 
        return final
    

# %%
