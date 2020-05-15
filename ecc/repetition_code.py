from encoder import *
from ecc_bio_interface import *
import random
def encode_messages(messages,poly_poses,degree):
    codewords = []
    for i in range(len(messages)):
        codewords.append(encode(degree,poly_poses,messages[i]))
    return codewords

def repeat_codewords(codewords,l):
    for i in range(l):
        codewords[i]+=codewords[i]
    return codewords

def join_for_synthesis(codewords):
    joined_codewords = [i for x in codewords for i in x]
    return joined_codewords

def map_codewords(codewords):
    return mapper(codewords)

def synthesise_and_shuffle(dna_seq_to_synthesise,n):
    shuffle = []
    for i in range(len(dna_seq_to_synthesise)):
        if i == len(dna_seq_to_synthesise)-n+1:
            break
        shuffle.append(dna_seq_to_synthesise[i:i+n])
    return shuffle

def remap_shuffle(shuffle):
    return [re_mapper(x) for x in shuffle]

def transpose_vector(shuffled_codewords):
    doubled_array = []
    for i in shuffled_codewords:
        outer_row=[]
        for j in i:
            inner_row=[]
            inner_row.append(j)
            outer_row.append(inner_row)
        doubled_array.append(outer_row)
    return doubled_array

def protonize(s,proto_messages):
	F4=[1,2,3]
	final_messages = []
	s_part = [1,2, 3,1,3]
	for message in proto_messages:
		s_part_new = []
		s_part_new =s_part + message
		final_messages.append(s_part_new)
		print(s_part)
		print(message)

	print("Final",final_messages)
	return final_messages
