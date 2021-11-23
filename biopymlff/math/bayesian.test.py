import numpy as np
import math

from matplotlib import pyplot

from biopymlff.math.bayesian import get_alpha_beta, get_sigma

x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
y = np.array([10, 20, 30, 40, 50, 60, 70 ,80, 90, 100])

# get_alpha_beta(x, y, 1, 1, 1000)
sigma = get_sigma(x, y)
print("sigma " + str(sigma))