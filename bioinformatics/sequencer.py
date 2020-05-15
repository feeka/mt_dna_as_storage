from interfaces.ecc_bio_interface import re_mapper,mapper
class Sequencer(object):
    """We are sequencing DNA """


    def create_reads(self,dna_seq_to_synthesise,n):
        """ Create ideal Genome reads out of given Genome
        :param dna_seq_to_synthesise: pass sequence to synthesize
        :param n: quantity of codewords
        """
        
        shuffle = []
        for i in range(len(dna_seq_to_synthesise)):
            if i == len(dna_seq_to_synthesise)-n+1:
                break
            shuffle.append(dna_seq_to_synthesise[i:i+n])
        self.reads = shuffle


    def remap_shuffle(self):
        self.shuffle = [re_mapper(x) for x in self.reads]

    def convert_to_string(self,valid_reads):
        self.valid_reads = valid_reads
        self.reads_of_string=[]
        for element in valid_reads:
            row=[]
            converted_to_string =""
            row=mapper(element)
            for i in row:
                converted_to_string +=i
            self.reads_of_string.append(converted_to_string)




    def create_k_mers(self,valid_reads,k):
        self.convert_to_string(valid_reads)
        k_mer=[]
        for read in self.reads_of_string:
            for j in range(len(read)-k+1):
                element = read[j:j+k]
                k_mer.append(element)
        self.k_mers =list(k_mer)
