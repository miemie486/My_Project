# x''+b*x'+A*x^3=B*cos(omega*t)
#A=1 omega=1 (B,b)=(7,6) 
# t=0 x=3 x'=0

import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
import multiprocessing
import time
from numba import jit
def df2(x,z,t,):
    a=B*m.cos(t)-x*x*x-b*z
    return a
kk=[[7,6],[7,0.6],[10,0.05],[7,0.01]]
p=0
while p<4:
    B=kk[p][0]
    b=kk[p][1]
    n=100000
    i=0
    h=25*m.pi/n
    r=np.array([3,0])
    xx=[]
    yy=[]
    while i<=n:
        k1=np.array([r[1],df2(r[0],r[1],i*h)])
        #print(df2(r[0],r[1],i*h))
        r2=r+0.5*k1*h
        k2=np.array([r2[1],df2(r2[0],r2[1],i*h+0.5*h)])
        r3=r+0.5*k2*h
        k3=np.array([r3[1],df2(r3[0],r3[1],i*h+0.5*h)])
        r4=r+k3*h
        k4=np.array([r4[1],df2(r4[0],r4[1],i*h+h)])
        rplus=r+(1/6)*(k1+2*k2+2*k3+k4)*h
        r=rplus
        xx.append(r[0])
        yy.append(i*h)
        i=i+1
    plt.figure(figsize=(10,10))
    plt.plot(yy,xx)
    a1='b_time'
    a2=str(kk[p][0])
    aa='_'
    a3=str(kk[p][1])
    a4='.png'
    plotPath= a1+a2+aa+a3+a4
    plt.savefig(plotPath)
    p=p+1
plt.show()
