from numba import jit
import numpy as np
c=12456
@jit
def f(x):

  a=np.array([1,2,3])
  a[2]=11
  for i in a:
    print(i)
  for i in x:
    print(i)
  print(c)
  return a


q=np.array([4,5,6])
print(f(q))
