# %%

def exon(a):
    
    dic ={
    "pc9_exon18" : "tcttgtcccccccagcttgtggagcctcttacacccagtggagaagctcccaaccaagctctcttgaggatcttgaaggaaactgaattcaaaaagatcaaagtgctgggctccggtgcgttcggcacggtgtataaggtaaggtccc".upper(),

    "pc9_exon19" : "ctctctctgtcatagggactctggatcccagaaggtgagaaagttaaaattcccgtcgctatcaaaacatctccgaaagccaacaaggaaatcctcgatgtgagtttct".upper(),

    #PC9 Exon19 deletion 15bp  
    # 109bp - 15bp ; 94bp

    "pc9_exon20" : "tctccctccctccaggaagcctacgtgatggccagcgtggacaacccccacgtgtgccgcctgctgggcatctgcctcacctccaccgtgcagctcatcacgcagctcatgcccttcggctgcctcctggactatgtccgggaacacaaagacaatattggctcccagtacctgctcaactggtgtgtgcagatcgcaaaggtaatcaggg".upper(),

    "pc9_exon21" : "tcttctctgtttcagggcatgaactacttggaggaccgtcgcttggtgcaccgcgacctggcagccaggaacgtactggtgaaaacaccgcagcatgtcaagatcacagattttgggctggccaaactgctgggtgcggaagagaaagaataccatgcagaaggaggcaaagtaaggaggt".upper(),

    "pc9_exon22" : "tctcaccatcccaaggtgcctatcaagtggatggcattggaatcaattttacacagaatctatacccaccagagtgatgtctggagctacggtgagtcata".upper(),

    "pc9_exon23" : "gcttcatcctctcaggggtgaccgtttgggagttgatgacctttggatccaagccatatgacggaatccctgccagcgagatctcctccatcctggagaaaggagaacgcctccctcagccacccatatgtaccatcgatgtctacatgatcatggtcaagtgtgagtgact".upper(),

    "pc9_exon24" : "tcattccttccccaggctggatgatagacgcagatagtcgcccaaagttccgtgagttgatcatcgaattctccaaaatggcccgagacccccagcgctaccttgtcattcaggtacaaattg".upper(),
    
    "pc9_exon2" : "cacaggtaagccaagatgggtgcatacaagtacatccaggagctatggagaaagaagcagtctgatgtcatgcgctttcttctgagggtccgctgctggcagtaccgccagctctctgctctccacagggctccccgccccacccggcctgataaagcgcgccgactgggctacaaggccaagcaaggtacgtgatc".upper(),
    
    "pc9_exon197519": "ctctctctgtcatagggactctggatcccagaaggtgagaaagttaaaattcccgtcgctatcaaggaattaagagaagcaacatctccgaaagccaacaaggaaatcctcgatgtgagtttct".upper(),
    "pc9_exon197520": "tctccctccctccaggaagcctacgtgatggccagcgtggacaacccccacgtgtgccgcctgctgggcatctgcctcacctccaccgtgcaActcatcaTgcagctcatgcccttcggctgcctcctggactatgtccgggaacacaaagacaatattggctcccagtacctgctcaactggtgtgtgcagatcgcaaaggtaatcaggg".upper(),
    "pc9_exon1975202" :"tctccctccctccaggaagcctacgtgatggccagcgtggacaacccccacgtgtgccgcctgctgggcatctgcctcacctccaccgtgcaGctcatcaTgcagctcatgcccttcggctgcctcctggactatgtccgggaacacaaagacaatattggctcccagtacctgctcaactggtgtgtgcagatcgcaaaggtaatcaggg".upper(),
    "pc9_exon197521": "tcttctctgtttcagggcatgaactacttggaggaccgtcgcttggtgcaccgcgacctggcagccaggaacgtactggtgaaaacaccgcagcatgtcaagatcacagattttgggcGggccaaactgctgggtgcggaagagaaagaataccatgcagaaggaggcaaagtaaggaggt".upper(),
    "pc9_exon79020": "TCTCCCTCCCTCCAGGAAGCCTACGTGATGGCCAGCGTGGACAACCCCCACGTGTGCCGCCTGCTGGGCATCTGCCTCACCTCCACCGTGCAGCTCATCATGCAGCTCATGCCCTTCGGCTGCCTCCTGGACTATGTCCGGGAACACAAAGACAATATTGGCTCCCAGTACCTGCTCAACTGGTGTGTGCAGATCGCAAAGGTAATCAGGG",
    "pc9_exon790201":"TCTCCCTCCCTCCAGGAAGCCTACGTGATGGCCAGCGTGGACAACCCCCACGTGTGCCGCCTGCTGGGCATCTGCCTCACCTCCACCGTGCAGCTCATTATGCAGCTCATGCCCTTCGGCTGCCTCCTGGACTATGTCCGGGAACACAAAGACAATATTGGCTCCCAGTACCTGCTCAACTGGTGTGTGCAGATCGCAAAGGTAATCAGGG"
    }

    return dic["pc9_exon"+str(a)]




## Genetic code

def gene_code(b):
    dic = {
        "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S","TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W","CTT":"L","CTC":"L","CTA":"L",
        "CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
        "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T","AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
        "GTT":"V","GTC":"V","GTA":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A","GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"
    }

    return dic[b]
