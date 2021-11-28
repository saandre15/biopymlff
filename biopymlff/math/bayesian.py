import numpy as np
import math

from matplotlib import pyplot

from scipy import optimize

def get_alpha_beta(X=None, Y=None, alpha=None, beta=None, iter_nums=None):
    """ 
    Optimizes the alpha and beta 
    """
    M = 2
    N = len(X)
    J = np.zeros(shape=(iter_nums, 1))
    fai = np.zeros(shape=(N, 4))

    for i in range(0, N):
        for j in range(0, 3):
            fai[i][j] = X[i] ^ (j)

    X = fai
    value = np.linalg.eigvals( beta * X @ np.transpose(X))
    m = len(value)

    gamma = 0
    # Calulates the Eigenspectrum
    for j in range(0, m):
        gamma = gamma + (value[j] / (value[j] + alpha))


    alphas = []
    betas = []
    iteration = []

    for i in range(0, iter_nums):
        
        # beta = 0.6 * 10 **26

        A = beta * np.transpose(X) @ X + alpha
    
        
        # Weight Value
        mN = beta * np.linalg.pinv(A) @ np.transpose(X) @ Y
        
        alpha = gamma / (mN @ np.transpose(mN))
        val = Y - (X @ mN)
        beta = ( N - gamma ) / (np.transpose(val) @ val)
        

        EmN = (1/2) * beta * val @ np.transpose(val) + ((1/2) * alpha  * mN @ np.transpose(mN))

        J[i] = ( (1/2) * M * math.log(alpha) + (1/2) * math.log(beta) - EmN - (1/2) * math.log(np.linalg.det(A)) - (1/2) * N * math.log(2 * math.pi))
        
        print("alpha " + str(alpha))
        print("beta "  + str(beta))
        print("J " + str(J[i]))
        alphas.append(alpha)
        betas.append(beta)
        iteration.append(i)

    maxi_at = -1
    maxi = -math.inf
    counter = 0
    for j in J:
        if j > maxi: 
            maxi = j
            maxi_at = counter
        counter+=1

    alpha_min = alphas[maxi_at]
    beta_min = betas[maxi_at]
    j_min=J[maxi_at]


    # def f(args: list):
    #     alpha = args[0]
    #     beta = args[1] 
        
    #     print("test")
    #     A = beta * np.transpose(X) @ X + alpha
        
    #     # Weight Value
    #     mN = beta * np.linalg.pinv(A) @ np.transpose(X) @ Y
    #     val = Y - (X @ mN)
    #     EmN = (1/2) * beta * val @ np.transpose(val) + ((1/2) * alpha  * mN @ np.transpose(mN))
    #     print("ALPHA " + str(alpha))
    #     print("BETA " + str(beta))
    #     J = ( (1/2) * M * math.log(alpha) + (1/2) * math.log(beta) - EmN - (1/2) * np.log(np.linalg.det(A)) - (1/2) * N * math.log(2 * math.pi))
    #     print("J " + str(J))
        
    #     return J


    # val = optimize.minimize(lambda args: -f(args), (alpha_min, beta_min), method='L-BFGS-B', bounds=((1e-10, math.inf), (1e-10, math.inf)))
    # print(val)


    
    fig = pyplot.figure()
    ax = pyplot.axes(projection='3d')
    ax.plot_wireframe(np.log(alphas), np.log(betas), J, color='green')
    ax.set_title('log (a) vs log (B) vs log P(Y|a, B)')
    pyplot.show()


    return alpha_min, beta_min

def get_sigma(X=None, Y=None):
    alpha, beta = get_alpha_beta(X, Y, 1, 1, 1000)

    N = len(X)
    fai = np.zeros(shape=(N, 4))

    for i in range(0, N):
        for j in range(0, 3):
            fai[i][j] = X[i] ^ (j)

    X = fai

    A = beta * np.transpose(X) @ X + alpha
    sigma = X[-1] @ np.linalg.pinv(A) @ np.transpose(X[-1])

    
    return sigma

def get_m_n(X=None, Y=None):
    raise NotImplementedError("Get M N for baysian statistics is not implemented.")
    
    
x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
y = np.array([10, 20, 30, 40, 50, 60, 70 ,80, 90, 100])

# get_alpha_beta(x, y, 1, 1, 1000)
sigma = get_alpha_beta(x, y, 1, 1, 400)