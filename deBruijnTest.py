from bioinformatics.dbg import *


# see if we can get correct reconstruction at k=4
# ACGGTCA
#  CGGTC
#   GGTCA

def split(word):
    return [char for char in word]

seq=['CCCGGC', 'CGAACT', 'TCCCCC', 'CCCCGT', 'GATAAC', 'ACTATT', 'ATACCA', 'AAGCCC', 'GCCCCC', 'CGTAGC', 'CCCTCC', 'CCCGAA', 'CCCAAA', 'CCAGGA', 'CCCACC', 'TGCAAA', 'CCACCC', 'AGCGAG', 'ATTCGC', 'CCCCAG', 'AGAACC', 'CTCCCC', 'TCGCAG', 'CCCCCA', 'ATAACG', 'GAACCC', 'CCCCCG', 'TAACGC', 'CACCCG', 'AAACCC', 'GAAGAC', 'CTATCC', 'ATCCCC', 'CCGAAG', 'CTGTGC', 'CCGGAT', 'TTATAC', 'CAGCCC', 'CCCCCC', 'AACCCC', 'CCCCGG', 'GTGCAA', 'CAATTC', 'CTAGAA', 'CCAGAA', 'CCGGCT', 'ACCCGG', 'ACTATC', 'CCCCAC', 'CGAAGC', 'CAAACC', 'CAGAAC', 'CCCAGG', 'CCCAGC', 'GCCGAA', 'TTCGCA', 'TACCAT', 'CCCTAG', 'CCCAAT', 'CTATTC', 'AGGATA', 'CCCTTA', 'CCCGTA', 'GCTTCC', 'CCTCCC', 'ACCATC', 'CCCCAA', 'CCCTGT', 'TATTCC', 'ACCCCC', 'CCCCGA', 'CCCGGA', 'GAAGCC', 'GACCCC', 'AATTCG', 'CGCCCC', 'CCGAAC', 'GCAAAC', 'CGCAGC', 'GGCTTC', 'AGACCC', 'TAGCGA', 'CCCAGA', 'AGAACT', 'CCAATT', 'TAGAAC', 'TATACC', 'CCAGCC', 'GAACTA', 'CAGGAT', 'CCGTAG', 'AGCCCC', 'CCTAGA', 'ACCCAC', 'AACTAT', 'CCCCCT', 'CCAAAC', 'AAACTA', 'TATCCC', 'CACCCC', 'TGTGCA', 'AAGACC', 'GGCCGA', 'AACGCC', 'GGATAA', 'CGAAGA', 'CCCCTA', 'CTTATA', 'TTCCCC', 'CTTCCC', 'CCCCTG', 'CCCCTT', 'CCATCC', 'CCCCTC', 'CAAACT', 'AACCCA', 'CCGGCC', 'GCAGCC', 'CATCCC', 'ACGCCC', 'CGGCTT', 'GTAGCG', 'CGGATA', 'CCTTAT', 'CCTGTG', 'CGGCCG', 'ATTCCC']
sequence=""
for i in sequence:
    sequence+=i
#norm = "ACGGTCGGAT"
a=set()
for i in range(len(sequence)-5):
    a.add(sequence[i:i+5])
seq=list(a)
print(seq)
walk = DeBruijnGraph(seq,5).eulerianWalkOrCycle()
res = walk[0] + ''.join(map(lambda x: x[-1], walk[1:]))
print(res)
