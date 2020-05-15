DNA_TO_ECC = {'A':0,'C':1,'G':2,'T':3}
def re_mapper(information):
	ecc_message = []
	for i in information:
		ecc_message.append(DNA_TO_ECC.get(i))
	return ecc_message


ECC_TO_DNA = {0:'A',1:'C',2:'G',3:'T'}
def mapper(information):
	dna_seq = []
	for i in information:
		dna_seq.append(ECC_TO_DNA.get(i))
	return dna_seq
