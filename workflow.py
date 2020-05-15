from ecc.encoder import Encoder
from bioinformatics.synthesizer import Synthesizer
from bioinformatics.sequencer import Sequencer
from ecc.codeword_detector import CodewordDetector
from bioinformatics.dbg import DeBruijnGraph, DeBruijnPlot
from collections import OrderedDict
from interfaces.ecc_bio_interface import *
from bioinformatics.scs import *
from interfaces.work_with_data import *

#######RAW DATA FOR TESTING########
proto_message1 = [3, 2, 1, 1, 1, 2, 1, 2]
proto_message2 = [1, 3, 2, 2, 1, 1, 2, 2]
proto_message3 = [1, 1, 2, 1, 3, 1, 3, 1]
proto_message4 = [2, 1, 1, 3, 1, 2, 2, 1]
proto_message5 = [2, 1, 2, 1, 2, 3, 1, 2]
proto_message6 = [1, 3, 1, 2, 1, 3, 1, 2]
k = len(proto_message1)
proto_messages = []
proto_messages.append(proto_message1)
proto_messages.append(proto_message2)
proto_messages.append(proto_message3)
proto_messages.append(proto_message4)
proto_messages.append(proto_message5)
proto_messages.append(proto_message6)
polynom = [1,1,1,2,3]
degree = 3
s=3
#######METHOD COMING SOON########



#--------------WE ENCODE MESSAGES-------------------#
ENCODER = Encoder(proto_messages,degree,polynom,s,k+s,15)
print("\n#------------------ENCODER----------------------#\n")
s_part=[1,2,3]
ENCODER.protonize(s_part)       #protonize given payload
ENCODER.encode_messages() #encode the resulting messages

print("\tResulting codewords:")
codewords = ENCODER.codewords
print("-----------------------------------------")

for i in range(len(codewords)):
    print("C"+str(i),codewords[i])

systematic_generator = ENCODER.systematic_generator

parity_check = ENCODER.parity_check
print("\n\tSystematic Generator Matrix:")
print("-----------------------------------------")
for i in systematic_generator:
    print(i)
print("\n#------------------ENCODER----------------------#\n")

#--------------WE ENCODE MESSAGES-------------------#

#--------------WE SYNTHESIZE CODEWORDS--------------#
print("\n#-------------SYNTHESIZER-------------------#")
SYNTHESIZER = Synthesizer(codewords)
print("\n\tSynthesized sequence:")
print("-----------------------------------------")
SYNTHESIZER.map_codewords()
print(SYNTHESIZER.mapped_to_dna)
mapped_to_dna = SYNTHESIZER.mapped_to_dna
print("\n#-------------SYNTHESIZER-------------------#\n")
#--------------WE SYNTHESIZE CODEWORDS--------------#

#----------------CREATE READS OUT OF CODEWORDS--------#
print("#----------------CREATE READS OUT OF CODEWORDS--------#")
SEQUENCER = Sequencer()

SEQUENCER.create_reads(mapped_to_dna,15)

pr = convert_to_string_reads(SEQUENCER.reads)
distort_reads(2,pr)
SEQUENCER.reads=convert_back_to_chars(pr)
SEQUENCER.remap_shuffle()

shuffle =SEQUENCER.shuffle
print("\n\tIDEAL Reads")
print("-----------------------------------------")
for i in SEQUENCER.shuffle: print(i)
print("#----------------CREATE READS OUT OF CODEWORDS--------#")
#----------------CREATE READS OUT OF CODEWORDS--------#


#----------------CODEWORD DETECTOR---------------------#
print("#----------------CODEWORD DETECTOR----------------#")
CD = CodewordDetector(parity_check,shuffle)
CD.perform_calculation_to_check()
print("\n\tPOSITIVES")
print("-----------------------------------------")
for i in CD.list_of_positives: print(i)
print("\n\tNEGATIVES")
print("-----------------------------------------")
for i in CD.list_of_negatives: print(i)
SEQUENCER.create_k_mers(CD.list_of_positives,4)
res = list(OrderedDict.fromkeys(SEQUENCER.k_mers))
previous_element=""
for i in range(len(res)):
    print(i*" ",res[i])
    previous_element=res[i][1:k-1]
original = ""
for i in SYNTHESIZER.mapped_to_dna:
	original+=i
print(len(original))
print(res)

count = 0
for i in SEQUENCER.reads_of_string:
    count+=1
    print(count*" ",i)

#############DBG test#####################


print(pr)
#walk=DeBruijnGraph(SEQUENCER.reads_of_string,2).eulerianWalkOrCycle()
from_distorted_reads = greedy_scs(pr,14)
#something = walk[0] + ''.join(map(lambda x: x[-1], walk[1:]))



del(SEQUENCER.reads_of_string[2])
del(SEQUENCER.reads_of_string[3])
from_cws = greedy_scs(SEQUENCER.reads_of_string,s)
print("Distorted reads:\t",from_distorted_reads)
print("From CWs:\t\t",from_cws)
print("Original:\t\t",original)
