from interfaces.ecc_helper import *


class CodewordDetector():
    """docstring for ."""

    def __init__(self, h_m,codewords):

        self.parity_check = h_m
        self.codewords=codewords

    def check_whether_codeword(self,list1):
    	for x in list1:
    		if x!= 0:
    			return False
    	return True

    def perform_calculation_to_check(self):
    	self.list_of_positives=[]
    	self.list_of_negatives=[]
    	count = 0
    	for vector in self.codewords:
            res=vec_mat(vector,transpose(self.parity_check))

            if self.check_whether_codeword(res):
                count+=1

                self.list_of_positives.append(vector)
            else:
                self.list_of_negatives.append(vector)
