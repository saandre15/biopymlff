import numpy as np
import math

from matplotlib import pyplot

def get_alpha_beta(X=None, Y=None, alpha=None, beta=None, iter_nums=None):
    M = 5
    N = len(X)
    J = np.zeros(shape=(iter_nums, 1))
    fai = np.zeros(shape=(N, 4))

    for i in range(0, N):
        for j in range(0, 3):
            fai[i][j] = X[i] ^ (j)

    X = fai
    value = np.linalg.eigvals( beta * X @ np.transpose(X))
    m = len(value)


    alphas = []
    betas = []
    iteration = []

    for i in range(0, iter_nums):
        gamma = 0
        # Calulates the Eigenspectrum
        for j in range(0, m):
            gamma = gamma + (value[j] / (value[j] + alpha))
        
        A = beta * np.transpose(X) @ X + alpha
    
        
        # Weight Value
        mN = beta * np.linalg.pinv(A) @ np.transpose(X) @ Y
        
        alpha = gamma / (mN @ np.transpose(mN))
        val = Y - (X @ mN)
        beta = ( N - gamma ) / (np.transpose(val) @ val)

        EmN = (1/2) * beta * val @ np.transpose(val) + ((1/2) * alpha  * mN @ np.transpose(mN))

        J[i] = ( (1/2) * M * math.log(alpha) + (1/2) * math.log(beta) - EmN - (1/2) * math.log(np.linalg.det(A)) - (1/2) * N * math.log(2 * math.pi))
        
        alphas.append(alpha)
        betas.append(beta)
        iteration.append(i)
    
    
    fourier = np.fft.fft(a=alphas[-21:-1])
    freq = np.fft.fftfreq(20, d=0.1)

    pyplot.plot(freq, [math.sqrt(math.pow(val.real, 2) + math.pow(val.imag, 2)) for val in fourier])
    pyplot.xlabel("Frequency")
    pyplot.ylabel("Amplitude")

    pyplot.show()

    pyplot.plot(iteration, alphas)
    pyplot.xlabel("Iterations")
    pyplot.ylabel("alpha(a)")
    pyplot.show()

    return alpha, beta

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
sigma = get_sigma(x, y)
print("sigma " + str(sigma))