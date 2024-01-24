from numba import jit
import numpy as np



@jit(nopython=True)
def f(Genes_arrary):
    fitness=0
    for i in Genes_arrary:
        fitness=fitness+i**2
    return fitness

Genes_arrary11=np.array([1,2,3])
f(Genes_arrary11)