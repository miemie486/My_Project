from cmath import inf
import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
import N_doping as N_D
import os
from mpl_toolkits.mplot3d import Axes3D
from timeit import default_timer as timer

from numba import cuda

@cuda.jit
def gpu():
    print(cuda.threadIdx.x)

@cuda.jit
def gpu_add(a,b,result,n):
    idx=0
    idx=cuda.threadIdx.x+cuda.blockDim.x*cuda.blockIdx.x
    if idx<n:
        for i in range(1000000):
            result[idx]=a[idx]+b[idx]


n=2000000
x=np.arange(n).astype(np.float32)
y=2*x
gpu_result=np.zeros(n)
cpu_result=np.zeros(n)

thread_per_block=1024
blocks_per_grid=int(m.ceil(n/thread_per_block))

start = timer()

gpu_add[blocks_per_grid,thread_per_block](x,y,gpu_result,n)

cuda.synchronize()
print("gpu",timer()-start)
start = timer()

cpu_result=np.add(x,y)

print("cpu",(timer()-start)*1000000)
