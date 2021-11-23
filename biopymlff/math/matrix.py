import numpy as np

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
    sum_sq = np.sum(matrix_sq)
	
    frob_col = np.sum(matrix_sq, axis=0)
    probs_col = frob_col/sum_sq				#probability of each column.
	
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
	
    X, sigma, Y_t = mysvd(W,rank)					#svd decomposition of W.
	
    for i in range(rank):
	    if(sigma[i,i] >= 1):						#taking pseudo-inverse of sigma.
		    sigma[i,i] = 1/sigma[i,i]
	    else:
		    sigma[i,i] = 0
	
    U = np.dot(Y_t.T, np.dot(np.dot(sigma,sigma), X.T))		#finding U.
	
    return C, U, R