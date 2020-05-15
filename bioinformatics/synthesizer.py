from interfaces.ecc_bio_interface import *


class Synthesizer(object):
    """Class serves as Synthesizer"""


    def __init__(self,codewords):
        self.codewords=codewords

    def __join_for_synthesis(self,codewords):
        """Concatenating to synthesize information
            :param codewords: pass codewords to synthesize
        """
        joined_codewords = [i for x in codewords for i in x]
        return joined_codewords

    def map_codewords(self):
        """Mapping to Genome to synthesize
            :param sequence: pass sequence to synthesize
        """
        self.joined_codewords=self.__join_for_synthesis(self.codewords)
        self.mapped_to_dna=mapper(self.joined_codewords)
