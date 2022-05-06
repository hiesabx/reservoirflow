from openresim import grids, fluids, models, wells
import scipy.sparse as ss
import numpy as np


# arr0 = np.ones((1000000,1000000), dtype='int')
arr = ss.lil_matrix([20,20], shape=(10,10), dtype='int') 
arr[0,0] = 10

print(arr.toarray())