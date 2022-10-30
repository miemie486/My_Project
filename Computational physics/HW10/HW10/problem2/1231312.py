import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
def f(x):
    a=(m.sin(m.sqrt(100*x)))**2
    return a
def f1(x):
    h=1/1000
    a=(f(x+h)-f(x-h))/(2*h)
    return a
def f2(x):
    h=1/1000
    a=(f1(x+h)-f1(x-h))/(2*h)
    return a
def f3(x):
    h=1/1000
    a=(f2(x+h)-f2(x-h))/(2*h)
    return a
def f4(x):
    h=1/1000
    a=(f3(x+h)-f3(x-h))/(2*h)
    return a
n=4
xx=np.linspace(n/1000,1,1001-n)
i=n
y=[]
while i<=1000:
    y.append(f4(i/1000))
    i=i+1
y=np.array(y)
plt.plot(xx,y)
plt.show()
