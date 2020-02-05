import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

k=10001
I=np.ones(k)
x=np.linspace(0,1,k)
y=924*x**6-2772*x**5+3150*x**4-1680*x**3+420*x**2-42*x+I
plt.plot(x,y)
plt.plot(x,0*x)
plt.show()