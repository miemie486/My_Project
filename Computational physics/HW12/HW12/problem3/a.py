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
    a=-y/2
    return a

def f(x):
    a=m.exp(-x/2)
    return a
@jit
def g():
    n=15
    while n<=15*10**2:
        h=15/n
        i=0
        y=1
        yy=[]
        xx=[]
        yyy=[]
        xxx=[]
        while i<=n:
            xx.append(i*h)
            xxx.append(i*h)
            yy.append(y)
            yyy.append(abs(f(i*h)-y))
            k1=df(i*h,y)
            yplus=y+k1*h
            #print(yplus)
            y=yplus
            i=i+1
        c1=str(h)
        c2='h='
        c3='error h='
        plt.figure(1)
        plt.plot(xx,yy,label=c2+c1)
        plt.legend()
        plt.figure(2)
        plt.plot(xxx,yyy,label=c3+c1)
        plt.legend()
        error=abs((f(1+h)-y)/f(1+h))
        #print(m.log10(n),error,f(1+h))
        n=n*10
        
g()
plt.show()
end=time.time()
print(end-s1)


