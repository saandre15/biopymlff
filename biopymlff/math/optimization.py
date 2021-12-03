import numpy as np
import os

from ase.atoms import Atoms
from ase.io.proteindatabank import read_proteindatabank

from ..descriptors.soap import SOAP

def svd_approximation(matrix, rank):
    
	m = matrix.shape[0]
	n = matrix.shape[1]
	
	if((rank>m) or (rank>n)):
		print("error: rank greater than matrix dimensions.\n")
		return;
		
	matrix_t = matrix.T
	
	A = np.dot(matrix, matrix_t)						#calculate matrix multiplied by its transpose
	values1, v1 = np.linalg.eigh(A)						#get eigenvalues and eigenvectors
	v1_t = v1.T
	v1_t[values1<0] = 0						#discarding negative eigenvalues and corresponding eigenvectors(they are anyway tending to zero)
	v1 = v1_t.T
	values1[values1<0] = 0
	#values1 = np.absolute(values1)
		
	values1 = np.sqrt(values1)						#finding singular values.
	
	idx = np.argsort(values1)						#sort eigenvalues and eigenvectors in decreasing order
	idx = idx[: :-1]
	values1 = values1[idx]
	v1 = v1[:, idx]
	
	U = v1
	
	A = np.dot(matrix_t, matrix)						#calculate matrix transpose multiplied by matrix.
	values2, v2 = np.linalg.eigh(A)						#get eigenvalues and eigenvectors
	#values2 = np.absolute(values2)
	v2_t = v2.T
	v2_t[values2<0] = 0						#discarding negative eigenvalues and corresponding eigenvectors(they are anyway tending to zero)
	v2 = v2_t.T
	values2[values2<0] = 0
	
	values2 = np.sqrt(values2)						#finding singular values.
	
	idx = np.argsort(values2)						#sort eigenvalues and eigenvectors in decreasing order.
	idx = idx[: :-1]
	values2 = values2[idx]
	v2 = v2[:, idx]
	
	V = v2
	V_t = V.T										#taking V transpose.
	
	sigma = np.zeros((m,n))
	
	if(m>n):										#setting the dimensions of sigma matrix.
		
		sigma[:n, :] = np.diag(values2)
			
	elif(n>m):
		sigma[:, :m] = np.diag(values1)
		
	else:
		sigma[:, :] = np.diag(values1)
			
	if(m > rank):									#slicing the matrices according to the rank value.
		U = U[:, :rank]
		sigma = sigma[:rank, :]
	
	if(n > rank):
		V_t = V_t[:rank, :]
		sigma = sigma[:, :rank]
	
	check = np.dot(matrix, V_t.T)					
	#case = np.divide(check, values2[:rank])
	
	s1 = np.sign(check)
	s2 = np.sign(U)
	c = (s1==s2)
	
	for i in range(U.shape[1]):						#choosing the correct signs of eigenvectors in the U matrix.
		if(c[0, i]==False):
			U[:, i] = U[:, i]*-1
	
	
	return U, sigma, V_t
    
def cur_approximation(matrix, rank):
    """
	INPUT: matrix: user-rating matrix, rank: desired rank.
	OUTPUT: returns C, U and R resulting from the CUR decomposition of the matrix.

    Note
    ----
    Mahoney, Michael W., and Petros Drineas. "CUR matrix decompositions for improved data analysis." Proceedings of the National Academy of Sciences 106.3 (2009): 697-702.
	"""
    m = matrix.shape[0]
    n = matrix.shape[1]
	
    if((rank>m) or (rank>n)):
        print("error: rank greater than matrix dimensions.\n")
        return;
		
    C = np.zeros((m, rank))
    R = np.zeros((rank, n))
	

    matrix_sq = matrix**2
    
    remove = np.zeros(shape=(len(matrix_sq), 1))
    sum_sq = np.sum(matrix_sq)
	
    frob_col = np.sum(matrix_sq, axis=0)
    probs_col = frob_col/sum_sq				#probability of each column.
	# Modify this porition
    count=0
    count1=0
    temp = 0
    idx = np.arange(n)						#array of column indexes.
    taken_c = []
    dup_c = []


	
    while(count<rank):
        i = np.random.choice(idx, p = probs_col)	#choosing column index based on probability.
        count1 = count1+1
        if(i not in taken_c):
            C[:, count] = matrix[:, i]/np.sqrt(rank*probs_col[i])	#taking column after dividing it with root of rank*probability.
            count = count+1
            taken_c.append(i)
            dup_c.append(1)
        else:										#discarding the duplicate column and increasing its count of duplicates.
            temp = taken_c.index(i)
            dup_c[temp] = dup_c[temp]+1
			
    C = np.multiply(C, np.sqrt(dup_c))				#multiply columns by root of number of duplicates.
			
    frob_row = np.sum(matrix_sq, axis=1)
    probs_row = frob_row/sum_sq					#probability of each row.
	
    count=0
    count1=0
    idx = np.arange(m)							#array of row indexes.
    taken_r = []
    dup_r = []
	
    while(count<rank):
        i = np.random.choice(idx, p = probs_row)			#choosing row index based on probability.
        count1 = count1+1
        if(i not in taken_r):
            R[count, :] = matrix[i, :]/np.sqrt(rank*probs_row[i])		#taking row after dividing it with root of rank*probability.
            count = count+1
            taken_r.append(i)
            dup_r.append(1)
        else:
            temp = taken_r.index(i)							#discarding the duplicate row and increasing its count of duplicates.
            dup_r[temp] = dup_r[temp]+1
		
    R = np.multiply(R.T, np.sqrt(dup_r))				#multiply rows by root of number of duplicates.
    R = R.T
	
    W = np.zeros((rank, rank))
	
    for i, I in enumerate(taken_r):
	    for j, J in enumerate(taken_c):				#forming the intersection matrix W.
		    W[i, j] = matrix[I, J]
	
    X, sigma, Y_t = svd_approximation(W,rank)					#svd decomposition of W.
	
    for i in range(rank):
	    if(sigma[i,i] >= 1):						#taking pseudo-inverse of sigma.
		    sigma[i,i] = 1/sigma[i,i]
	    else:
		    sigma[i,i] = 0
	
    U = np.dot(Y_t.T, np.dot(np.dot(sigma,sigma), X.T))		#finding U.
	
    return C, U, R

def cur_approximation_modified(matrix, rank):
    """
	INPUT: matrix: user-rating matrix, rank: desired rank.
	OUTPUT: returns C, U and R resulting from the CUR decomposition of the matrix.

    Note
    ----
    Mahoney, Michael W., and Petros Drineas. "CUR matrix decompositions for improved data analysis." Proceedings of the National Academy of Sciences 106.3 (2009): 697-702.
	"""
    m = matrix.shape[0]
    n = matrix.shape[1]
	
    if((rank>m) or (rank>n)):
        print("error: rank greater than matrix dimensions.\n")
        return;
		
    R = np.zeros((rank, n))
	
    theshold = 10 ** -10

    matrix_sq = matrix**2

    
    eigenvalues = np.linalg.eig(matrix)[0]

    counter = 0
    for i in range(0, len(matrix_sq)):
        for j in range(0, len(matrix_sq[0])):
            eigenvalue = eigenvalues[j]
            if eigenvalue < theshold: matrix_sq[i][j] = 0
            
    sum_sq = np.sum(matrix_sq)
	
    frob_col = np.sum(matrix_sq, axis=0)
    print(frob_col)
    
    probs_col = frob_col/sum_sq				#probability of each column.
    print(probs_col)
    
    count=0
    count1=0
    temp = 0
    idx = np.arange(n)						#array of column indexes.
    taken_c = []
    dup_c = []
	
    rank_col = len(probs_col[probs_col != 0])
    C = np.zeros((m, rank_col))

    while(count < rank_col):
		# NOTE: The count will never be greater than the rank if the there are lenghts of the colomns that are not 0
        # print(idx)
        i = np.random.choice(idx, p = probs_col)	#choosing column index based on probability.
        # print(i)
        count1 = count1+1
        if(i not in taken_c):
            C[:, count] = matrix[:, i]/np.sqrt(rank*probs_col[i])	#taking column after dividing it with root of rank*probability.
            count = count+1
            taken_c.append(i)
            dup_c.append(1)
        else:										#discarding the duplicate column and increasing its count of duplicates.
            temp = taken_c.index(i)
            dup_c[temp] = dup_c[temp]+1
			
    C = np.multiply(C, np.sqrt(dup_c))				#multiply columns by root of number of duplicates.
			
    frob_row = np.sum(matrix_sq, axis=1)
    probs_row = frob_row/sum_sq					#probability of each row.
	
    count=0
    count1=0
    idx = np.arange(m)							#array of row indexes.
    taken_r = []
    dup_r = []
	
    while(count<rank):
	    i = np.random.choice(idx, p = probs_row)			#choosing row index based on probability.
	    count1 = count1+1
	    if(i not in taken_r):
		    R[count, :] = matrix[i, :]/np.sqrt(rank*probs_row[i])		#taking row after dividing it with root of rank*probability.
		    count = count+1
		    taken_r.append(i)
		    dup_r.append(1)
	    else:
		    temp = taken_r.index(i)							#discarding the duplicate row and increasing its count of duplicates.
		    dup_r[temp] = dup_r[temp]+1
		
    R = np.multiply(R.T, np.sqrt(dup_r))				#multiply rows by root of number of duplicates.
    R = R.T
	
    W = np.zeros((rank, rank))
	
    for i, I in enumerate(taken_r):
	    for j, J in enumerate(taken_c):				#forming the intersection matrix W.
		    W[i, j] = matrix[I, J]
	
    X, sigma, Y_t = svd_approximation(W,rank)					#svd decomposition of W.
	
    for i in range(rank):
	    if(sigma[i,i] >= 1):						#taking pseudo-inverse of sigma.
		    sigma[i,i] = 1/sigma[i,i]
	    else:
		    sigma[i,i] = 0
	
    U = np.dot(Y_t.T, np.dot(np.dot(sigma,sigma), X.T))		#finding U.
    return C, U, R

def oblique_projection_matrix(matrix: np.ndarray):
	dimension = matrix.shape
	rank = np.linalg.matrix_rank(matrix)
	U, S, V = svd_approximation(matrix, rank)
	B = matrix+np.identity(dimension[0])- (U @ U.T)
	return B

mol = read_proteindatabank(os.getcwd() + "/data/systems/4znn.not_wat.pdb")
# mol = Atoms(symbols='H2O')
fingerprint = SOAP(l_max=6, n_max=3, atom_sigma=0.4, r_cutoff=3, radial_scaling=1, cutoff_trans_width=1, central_weight=1,
    n_sparse=1, delta=1, covariance_type="dot_product", zeta=2.5, species=mol.get_chemical_symbols(), sparse_method="cur_points").to_tensor(mol)

print(cur_approximation(fingerprint, np.linalg.matrix_rank(fingerprint)))