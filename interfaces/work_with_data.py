from interfaces.ecc_bio_interface import mapper
import random

def distort_as_soa(dna,p):
    alphabet = ['C','T','G','A']
    for i in range(len(dna)):
        s=random.uniform(0,1)
        if s<p:
            dna[i] = random.choice(alphabet)



def convert_to_string_reads(mapped_to_dna):
    reads_string = []

    for reads in mapped_to_dna:
        read_string = ""
        for read in reads:
            read_string+=read
        reads_string.append(read_string)
    return reads_string

def distort_reads_strong_version(reads,p=0):
    t_reads = list(reads)
    everything_except_zero = ['C','T','G']
    everything_except_one = ['A','T','G']
    everything_except_two = ['C','A','T']
    everything_except_three = ['A','C','G']
    count = []
    cnt=0
    for i in range(len(t_reads)):
        if i>0:
            cnt+=1
        for j in range(len(t_reads[i])):

            s=random.uniform(0,1)
            if s<p:
                st = "["+str(cnt)+"]["+str(j)+"]"+t_reads[i][j]
                if t_reads[i][j]=='A':
                    del(t_reads[i][j])
                    t_reads[i].insert(j,random.choice(everything_except_zero))
                elif t_reads[i][j]=='C':
                    del(t_reads[i][j])
                    t_reads[i].insert(j,random.choice(everything_except_one))
                elif t_reads[i][j]=='G':
                    del(t_reads[i][j])
                    t_reads[i].insert(j,random.choice(everything_except_two))
                elif t_reads[i][j]=='T':
                    del(t_reads[i][j])
                    t_reads[i].insert(j,random.choice(everything_except_three))

                st=st+"->"+t_reads[i][j]
                count.append(st)


    #print("Positions distorted",count)
    return t_reads

def distort_reads(reads,p=0):
    t_reads = list(reads)
    inner_len =  len(reads[0])-1
    everything_except_zero = ['C','T','G']
    everything_except_one = ['A','T','G']
    everything_except_two = ['C','A','T']
    everything_except_three = ['A','C','G']
    count = []
    cnt=0
    for i in range(0,len(t_reads)):
        if i>0:
            cnt+=1
        j=random.randint(0,inner_len-1)
        s=random.uniform(0,1)
        if s<p:
            st = "["+str(cnt)+"]["+str(j)+"]"+t_reads[i][j]
            if t_reads[i][j]=='A':
                del(t_reads[i][j])
                t_reads[i].insert(j,random.choice(everything_except_zero))
            elif t_reads[i][j]=='C':
                del(t_reads[i][j])
                t_reads[i].insert(j,random.choice(everything_except_one))
            elif t_reads[i][j]=='G':
                del(t_reads[i][j])
                t_reads[i].insert(j,random.choice(everything_except_two))
            elif t_reads[i][j]=='T':
                del(t_reads[i][j])
                t_reads[i].insert(j,random.choice(everything_except_three))
            st=st+"->"+t_reads[i][j]
            count.append(st)

    #print("Positions distorted",count)
    return t_reads


def convert_back_to_chars(data):
    temp=[]
    for i in data:
        row=[]
        for j in i:
            row.append(j)
        temp.append(row)
    return temp

def generate_random_messages(num_of_chars):
    messages = []
    ffour = [0,1,2,3]
    for i in range(num_of_chars):
        messages.append(random.choice(ffour))
    return messages
