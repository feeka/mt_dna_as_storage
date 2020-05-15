from interfaces.helper_math import LookUpTable

class F_Four:
	elements = {0:[0,0],
				1: [0,1],
				2: [1,0],
				3: [1,1]}

	def __init__(self,n):
		self.n = n
		self.poly = self.elements.get(n)
		self.lut=LookUpTable()

	def perform_calculation(self,other,SIGN):
		if isinstance(other,int):
			other=F_Four(other)
		str_res = str(self.n)+ SIGN + str(other.n)

		result = self.lut.get_table(SIGN).get(str_res)
		return result

	def __add__(self,other):
		return self.perform_calculation(other,"+")

	def __sub__(self,other):
		return self.perform_calculation(other,"-")

	def __mul__(self,other):
		return self.perform_calculation(other,"*")

	def __truediv__(self,other):
		return self.perform_calculation(other,"/")

	def __str__(self):
		return str(self.n)

	def __eq__(self,other):
		if isinstance(other,int):
			other = F_Four(other)
		return self.n == other.n
