from interfaces.f_four import F_Four
from interfaces.ecc_bio_interface import *
"""Codeword detection

This script allows the user to detect whether code or not by applying various
algebraic operations

    * matrixmult - multiply 2 matrixes and return result in F_Four
    * transpose - transpose matrix and return result in F_Four
    * vec_mat - multiply vector by matrix  and return result in F_Four
    * rref - bring matrix to reduced row echelon form (helper function) - pass-by-value principle
    *
"""

def make_one_vector(multivector):
	one_vector = []
	for double_vector in multivector:
		element=0
		for j in double_vector:
			element=j
		one_vector.append(element)
	return one_vector


def matrixmult (A, B):
	rows_A = len(A)
	cols_A = len(A[0])
	rows_B = len(B)
	cols_B = len(B[0])
	if cols_A != rows_B:
		print("Cannot multiply the two matrices. Incorrect dimensions.")
		return
    # Create the result matrix
    # Dimensions would be rows_A x cols_B
	C = [[F_Four(0) for row in range(cols_B)] for col in range(rows_A)]
	A_new = []
	for vec in A:
		row = []
		for i in vec:
			row.append(F_Four(i))
		A_new.append(row)
	B_new = []
	for vec in B:
		row = []
		for i in vec:
			row.append(F_Four(i))
		B_new.append(row)

	for i in range(rows_A):
		for j in range(cols_B):
			#if isinstance(C[i][j],int):
			for k in range(cols_A):
				intermediate = F_Four(A_new[i][k] * B_new[k][j])
				C[i][j] = intermediate +C[i][j]
	return C

# transpose :: Matrix a -> Matrix a
def transpose(m):
    if m:
        inner = type(m[0])
        z = zip(*m)
        return (type(m))(
            map(inner, z) if tuple != inner else z
        )
    else:
        return m

def vec_mat(vec,M):
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


def rref( M):
    if not M: return
    lead = 0
    rowCount = len(M)
    columnCount = len(M[0])
    for r in range(rowCount):
        if lead >= columnCount:
            return
        i = r
        while M[i][lead] == F_Four(0):
            i += 1
            if i == rowCount:
                i = r
                lead += 1
                if columnCount == lead:
                    return

        M[i],M[r] = M[r],M[i]
        lv = M[r][lead]
        if isinstance(lv,int):
            lv = F_Four(lv)
        row = []
        for mrx in M[r]:
            if isinstance(mrx,int):
                mrx = F_Four(mrx)
            row.append(mrx/lv)
        M[r]=row
        for i in range(rowCount):
            if i != r:
                lv = M[i][lead]
                row1 = []
                if isinstance(lv,int):
                    lv=F_Four(lv)
                for rv, iv in zip(M[r],M[i]):
                    if isinstance(rv,int):
                        rv=F_Four(rv)
                    if isinstance(iv,int):
                        iv=F_Four(iv)
                    row1.append(iv - lv*rv)
                M[i] = row1

        lead += 1



def g_and_h(matrix,n,k):
    rref(matrix)
    matrix_A = []
    for i in range(len(matrix)):
        row =[]
        for j in range(k,n):
            row.append(matrix[i][j])
        matrix_A.append(row)
    matrix_H = transpose(matrix_A)
    counter = k
    for i in range(0,len(matrix_H)):
        row = []
        for j in range(k,n):
            if(counter==j):
                matrix_H[i].append(1)
                continue
            matrix_H[i].append(0)
        counter+=1

    return (matrix_H,matrix)
