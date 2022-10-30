import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
from numba import jit

#   y'=1-x+4y  y(0)=1
#   x from 0 to 10
s1=time.time()

def df(x,y):
    a=1-x+4*y
    return a

def f(x):
    a=x/4-3/16+(19/16)*m.exp(4*x)
    return a
@jit
def g():
    n=10
    while n<=10**2:
        h=1/n
        i=0
        y=1
        #yy=[]
        #xx=[]
        while i<=n:
            #xx.append(i*h)
            #yy.append(y)
            k1=df(i*h,y)
            k2=df(i*h+(3/4)*h,y+(3/4)*k1*h)
            yplus=y+((1/3)*k1+(2/3)*k2)*h
            #print(yplus)
            y=yplus
            i=i+1
        #plt.plot(xx,yy)
        #plt.show()
        error=abs((f(1+h)-y)/f(1+h))
        print(m.log10(n),error,f(1+h))
        n=n*10
        
g()
end=time.time()
print(end-s1)


