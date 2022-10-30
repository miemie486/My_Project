import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
from numba import jit

#   y'=1-x+4y  y(0)=1
#   x from 0 to 10
s1=time.time()

def df(x,y,z):
    a=-y/z
    return a
def df1(x,y,z,x1,y1,z1):
    a=y1/z1-y/z
    return a
def f(x,z):
    a=m.exp(-x/z)
    return a
@jit
def g():
    n=15*100
    for bkb in [0.02,2,200]:
        h=15/n
        i=0
        y=1
        y1=1
        yy=[]
        yy1=[]
        xx=[]
        xx1=[]
        yyy=[]
        xxx=[]
        while i<=n:
            xx.append(i*h)
            xx1.append(i*h)
            xxx.append(i*h)
            yy.append(y)
            yy1.append(y1)
            yyy.append(abs(f(i*h,bkb)-y))
            k1=df(i*h,y,2)
            yplus=y+k1*h
            k11=df1(i*h,y1,bkb,i*h,y,2)
            yplus1=y1+k11*h
            #print(yplus)
            y1=yplus1
            y=yplus
            i=i+1
        c1=str(bkb)
        c2='tau='
        c3='error tau='
        plt.figure(1)
        plt.plot(xx,yy,label=c2+c1)
        plt.figure(2)
        plt.plot(xxx,yyy,label=c2+c1)
        plt.legend()
        error=abs((f(1+h,bkb)-y)/f(1+h,bkb))
        #print(m.log10(n),error,f(1+h))
        
g()
plt.show()
end=time.time()
print(end-s1)


