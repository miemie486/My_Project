from numba import jit
import numpy as np
c=12456
@jit
def f(x):

  kk=[c,c,c]
  return kk


q=np.array([4,5,6])
print(f(q))
