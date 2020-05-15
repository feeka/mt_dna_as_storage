from ecc.encoder import Encoder
from bioinformatics.synthesizer import Synthesizer
from bioinformatics.sequencer import Sequencer
from ecc.codeword_detector import CodewordDetector
from bioinformatics.dbg import DeBruijnGraph, DeBruijnPlot
from collections import OrderedDict
from interfaces.ecc_bio_interface import *
from bioinformatics.scs import *
from interfaces.work_with_data import *


class WorkflowManager(object):
    """docstring for WorkflowManager
    """
    def __init__(self,proto_messages,polynom,s_part,degree,n):
        self.proto_messages = proto_messages
        self.degree = degree
        self.polynom = polynom
        self.s = len(s_part)
        self.l = len(proto_messages[0])
        self.n = n
        self.k = self.l+self.s
        self.s_part = s_part
        self.ENCODER = Encoder(self.proto_messages,self.degree,self.polynom,self.s,self.k,self.n)

    def generate_messages(self):
        #print("DATASET TEST STARTED")
        self.ENCODER.protonize(self.s_part)
        """s=input("Display PROTO MESSAGES or continue?")
        if s=="DR":
            print("\tProto Messages:")
            for i in range(len(self.ENCODER.messages)):
                print("message",i,self.ENCODER.messages[i])

        print("")
        print("GOING TO ENCODE MESSAGE")"""


    def encode(self):
        self.ENCODER.encode_messages() #encode the resulting messages
        self.parity_check = self.ENCODER.parity_check
        """s=input("Display CODEWORDS or continue?")
        if s=="DR":
            systematic_generator = self.ENCODER.systematic_generator
            print("\tSystematic Generator Matrix:")
            print("-----------------------------------------")
            for i in systematic_generator:
                print(i)
            print("")
            print("\tResulting codewords:")
            print("-----------------------------------------")

            for i in range(len(self.ENCODER.codewords)):
                print("C"+str(i),self.ENCODER.codewords[i])
        """
        print("")
        print("GOING TO MAP AND SYNTHESIZE INFORMATION")

    def map_and_synthesize(self):
        #--------------WE SYNTHESIZE CODEWORDS--------------#
        #s=input("Display CODEWORDS or continue?")
        self.SYNTHESIZER = Synthesizer(self.ENCODER.codewords)
        self.SYNTHESIZER.map_codewords()
        """if s=="DR":

            print("\n\tSynthesized DNA sequence")
            print(self.SYNTHESIZER.mapped_to_dna)

        print("")
        print("GOING TO SHUFFLE AND GENERATE READS FROM INFORMATION")
        """

    def create_reads_sequence(self,probability = 0):
        #----------------CREATE READS OUT OF CODEWORDS--------#
        #print("#----------------CREATE READS OUT OF CODEWORDS--------#")
        self.SEQUENCER = Sequencer()
        self.SEQUENCER.create_reads(self.SYNTHESIZER.mapped_to_dna,self.n)
        self.SEQUENCER.reads = distort_reads(self.SEQUENCER.reads,probability)
        self.probability_a=probability
        self.pr = convert_to_string_reads(self.SEQUENCER.reads)
        self.SEQUENCER.remap_shuffle()
        #s=input("Display READS or continue?")
        self.shuffle =self.SEQUENCER.shuffle
        """if s=="DR":
            print("\nReads with",distorted_position,"errors")
            print("-----------------------------------------")
            for i in self.SEQUENCER.shuffle: print(i)

        print("")
        print("GOING TO DETECT THE CODEWORDS")
        """
        #----------------CREATE READS OUT OF CODEWORDS--------#

    def detect_codewords(self):
        self.CD = CodewordDetector(self.parity_check,self.shuffle)
        self.CD.perform_calculation_to_check()
        """s=input("Display READS or continue?")
        if s=="DR":
            print("\n\tPOSITIVES")
            print("-----------------------------------------")
            for i in self.CD.list_of_positives: print(i)
            print("\n\tNEGATIVES")
            print("-----------------------------------------")
            for i in self.CD.list_of_negatives: print(i)
        print("")
        print("GOING TO ASSEMBLY PROBLEM")"""

    def assemble(self):
        self.SEQUENCER.create_k_mers(self.CD.list_of_positives,self.s)
        self.pr=convert_to_string_reads(self.SEQUENCER.reads)
        t=""
        for i in self.pr:
            t+=i
        #print(t)
        a=set()
        v=self.s-1
        s_mers=[]
        for i in range(len(t)-self.s):

            a.add(t[i:i+self.s])

        seq=list(a)
        #print(self.SEQUENCER.reads_of_string)
        #walk = DeBruijnGraph(seq,v).eulerianWalkOrCycle()
        #self.from_cws = walk[0] + ''.join(map(lambda x: x[-1], walk[1:]))
        #print(self.from_cws)
        #print("Reads",convert_to_string_reads(self.SEQUENCER.reads))

        self.from_cws = greedy_scs(self.SEQUENCER.reads_of_string,v)
        self.original = ""
        for i in self.SYNTHESIZER.mapped_to_dna:
        	self.original+=i

        """
        for i in range(len(self.from_cws),self.n):
            if self.s_part==self.from_cws[i+self.s]:
                print("equal")
                continue
            else:
                print("not")
                self.from_cws.replace(self.from_cws[i:i+self.s],self.s_part)
        """
        #print("Original",self.original)
        #print("From CWs",self.from_cws)
        return self.from_cws==self.original

    def compare_data(self):
        #self.from_distorted_reads = greedy_scs(self.pr,self.s)
        diff_num = abs(len(self.from_cws)-len(self.original))
        #self.from_cws=self.from_cws[:len(self.from_cws)-diff_num]
        u = zip(self.original,self.from_cws)
        #print("FROM CODEWORDS")
        #print("From dist reads",self.from_distorted_reads)
        count_for_cw = 0
        for i in range(len(self.from_cws)):
            if i==len(self.original):
                break
        #print(count_for_cw)
        count_for_distr = 0
        #print("DATASET TEST FINISHED")
        return (int(count_for_cw), 0)
