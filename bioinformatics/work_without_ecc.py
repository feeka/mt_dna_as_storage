from interfaces.ecc_bio_interface import *
import random
from bioinformatics.scs import *

class WorkWithoutECC:

    def __init__(self,pms,p):
        self.proto_messages = pms
        self.p = p

    def concatenate_to_one_string(self):
        strig = ""
        converted_to_dna = []
        for information in self.proto_messages:
            converted_to_dna.append(mapper(information))

        for i in converted_to_dna:
            for j in i:
                strig+=str(j)
        self.one_str = strig

    def distort(self):
        strig = ""
        everything_except_zero = ['C','T','G']
        everything_except_one = ['A','T','G']
        everything_except_two = ['C','A','T']
        everything_except_three = ['A','C','G']
        for i in self.one_str:
            s=random.uniform(0,1)
            if s<self.p:
                a=''
                if i=='A':
                    a=random.choice(everything_except_zero)
                elif i=='C':
                    a=random.choice(everything_except_one)
                elif i=='G':
                    a=random.choice(everything_except_two)
                elif i=='T':
                    a=random.choice(everything_except_three)
                strig+=a
            else:
                strig+=i
        self.distorted_str=strig


    def make_reads(self):
        k = len(self.proto_messages[0])
        reads = set()
        for i in range(len(self.one_str)):
            read_k=""
            if i<=len(self.distorted_str)-k:
                read_k+=self.distorted_str[i:i+k]
            else:
                break
            reads.add(read_k)
        self.reads=list(reads)

    def distort_reads_weak(self):
        everything_except_zero = ['C','T','G']
        everything_except_one = ['A','T','G']
        everything_except_two = ['C','A','T']
        everything_except_three = ['A','C','G']
        count = 0
        for i in range(len(self.reads)):
            s=random.uniform(0,1)
            #print("random",s)
            j=random.randint(0,len(self.reads[i])-1)
            if s<self.p:
                count+=1
                if self.reads[i][j]=='A':
                    self.reads[i] = self.reads[i][:j] + random.choice(everything_except_zero) + self.reads[i][j+1:]
                elif self.reads[i][j]=='C':
                    self.reads[i] = self.reads[i][:j] + random.choice(everything_except_one) + self.reads[i][j+1:]
                elif self.reads[i][j]=='G':
                    self.reads[i] = self.reads[i][:j] + random.choice(everything_except_two) + self.reads[i][j+1:]
                elif self.reads[i][j]=='T':
                    self.reads[i] = self.reads[i][:j] + random.choice(everything_except_three) + self.reads[i][j+1:]


    def distort_reads(self):
        everything_except_zero = ['C','T','G']
        everything_except_one = ['A','T','G']
        everything_except_two = ['C','A','T']
        everything_except_three = ['A','C','G']
        count=0
        for i in range(len(self.reads)):
            s=random.uniform(0,1)
            for j in range(len(self.reads[i])):
                #s=random.uniform(0,1)

                if s<self.p:

                    count+=1
                    if self.reads[i][j]=='A':
                        self.reads[i] = self.reads[i][:j] + random.choice(everything_except_zero) + self.reads[i][j+1:]
                    elif self.reads[i][j]=='C':
                        self.reads[i] = self.reads[i][:j] + random.choice(everything_except_one) + self.reads[i][j+1:]
                    elif self.reads[i][j]=='G':
                        self.reads[i] = self.reads[i][:j] + random.choice(everything_except_two) + self.reads[i][j+1:]
                    elif self.reads[i][j]=='T':
                        self.reads[i] = self.reads[i][:j] + random.choice(everything_except_three) + self.reads[i][j+1:]

    def compare_data(self):

        from_dist_reads = greedy_scs(list(self.reads),len(self.reads[0])-1)
        u=zip(self.one_str,from_dist_reads)
        count_for_distr = 0
        for i,j in u:
            if i!=j:
                count_for_distr +=0.8
        #diff_num = abs(len(from_dist_reads)-len(self.one_str))

        return count_for_distr
