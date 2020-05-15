from interfaces.f_four import F_Four
from interfaces.ecc_helper import g_and_h

class Encoder(object):
	"""Class serves as an Encoder"""
	def __init__(self,proto_messages,degree,polynom,s,k,n):
		self.proto_messages = proto_messages
		self.degree=degree
		self.polynom=polynom
		self.s=s
		self.k=k
		self.n=n


	def protonize(self,s_part):
		F4=[1,2,3]
		final_messages = []

		for message in self.proto_messages:
			s_part_new = []
			s_part_new =s_part + message
			final_messages.append(s_part_new)
		self.messages =  final_messages

	def vec_mat(self,vec,M):
		vec_res =[]
		transposed =list(zip(*M))
		transposed_new = []
		for vecs in transposed:
			els = []
			for i in vecs:
				els.append(F_Four(i))
			transposed_new.append(els)
		vec_new = []
		for i in vec:
			vec_new.append(F_Four(i))
		for i in range(len(transposed_new)):
			summa=F_Four(0)
			for j in range(len(transposed_new[0])):
				umnoj = vec_new[j] * transposed_new[i][j]
				summa = summa + F_Four(umnoj)
				summa = F_Four(summa)
			vec_res.append(summa.n)
		return vec_res


	def __construct_generator_matrix(self,msg_length,polynom):
		generator_matrix = []
		for i in range(msg_length):
			generator_matrix.append(polynom)
		g_m =[]
		for i in range(len(generator_matrix)):
			g_m.append([0]*i + generator_matrix[i] + [0]*(msg_length-i-1))
		return g_m

	def encode(self,message):
		g_matr = self.__construct_generator_matrix(self.k,self.polynom)
		g2=g_and_h(g_matr,self.n,self.k)
		self.systematic_generator = g2[1]
		code = self.vec_mat(message,g2[1])
		self.parity_check = g2[0]
		return code

	def encode_messages(self):
		codewords = []
		for i in range(len(self.messages)):
			codewords.append(self.encode(self.messages[i]))
		self.codewords=codewords
