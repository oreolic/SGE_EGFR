from itertools import combinations


def _change_mm_num_to_nt(i,j,mm_bc_dic,bc,mmdic):
    wt_nt = bc[j]
    mm_nt_lst = mmdic[wt_nt]

    new_dic = {}
    for mmbcs in mm_bc_dic:
        for mm_nt in mm_nt_lst:
            ed_bc = mmbcs.replace(str(i),mm_nt)
            new_dic[ed_bc] = bc

    new_dic
    return new_dic


def _generate_mm_bc(sb,n):
    poslst = [i for i in range(len(sb))]
    mm_comb = list(combinations(poslst,n))
    mmdic = {'A':['T','G','C'],'T':['A','G','C'],'G':['A','T','C'],'C':['A','T','G']}       

    final_dic = {sb:sb}

    for tup in mm_comb:
        mmbc = ''

        m = 0 
        for idx,nt in enumerate(sb): ## mmbc : TGT0TCTTTC1CGTATAC
            if idx in tup:
                mmbc += str(m)
                m += 1
            else:
                mmbc += nt

        mm_bc_dic = {mmbc:1}
        for i,j in enumerate(tup):
            mm_bc_dic = _change_mm_num_to_nt(i,j,mm_bc_dic,sb,mmdic)

        for mm_ed_bc in mm_bc_dic:
            final_dic[mm_ed_bc] = sb
    return final_dic


def n_mm_barcode(sb,n):
    if n == 0:
        merged_dic = {sb:sb}
    else:
        merged_dic = {}

        for i in range(1,n+1):
            eachdic = _generate_mm_bc(sb,i)
            for ed_bc in eachdic:
                merged_dic[ed_bc] = sb
    return merged_dic


