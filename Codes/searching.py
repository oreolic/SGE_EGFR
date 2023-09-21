

def search_motif(read,searching_length,ed_dic,start=0,end = True):
    n= start
    motif = 'FALSE'
    motif_index = -1
    if end == True:
        end = len(read)-searching_length+1
    else:
        end = end

    while n < end:
        seq = read[n:n+searching_length]
        if seq in ed_dic:
            motif = ed_dic[seq]
            motif_index = n
            break
        else: 
            n += 1

    return motif, motif_index


