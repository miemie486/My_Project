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
s1=time.time()
B=7
b=0.01

def g(oo):
    n=1000000
    h=10000*m.pi/n
    r=np.array([3,0])
    xx=[]
    yy=[]
    for i in range(0,n):
        k1=np.array([r[1],df2(r[0],r[1],i*h,oo)])
        #print(df2(r[0],r[1],i*h))
        r2=r+0.5*k1*h
        k2=np.array([r2[1],df2(r2[0],r2[1],i*h+0.5*h,oo)])
        r3=r+0.5*k2*h
        k3=np.array([r3[1],df2(r3[0],r3[1],i*h+0.5*h,oo)])
        r4=r+k3*h
        k4=np.array([r4[1],df2(r4[0],r4[1],i*h+h,oo)])
        rplus=r+(1/6)*(k1+2*k2+2*k3+k4)*h
        r=rplus
        if(i%200==0):
            xx.append(r[0])
            yy.append(r[1])
    plt.figure(figsize=(10,10))
    plt.scatter(xx,yy)
    a1='c_phase'
    a2=str(oo)
    a3='.png'
    plotPath=a1+a2+a3
    plt.savefig(plotPath)

def df2(x,z,t,oo):
    a=B*m.cos(t+oo*(m.pi/360))-x*x*x-b*z
    return a

def func1(we):
    for oo in range(0,91):
        g(oo)
    print(1)
    
def func2(we):
    for oo in range(91,181):
        g(oo)
    print(2)

def func3(we):
    for oo in range(181,271):
        g(oo)
    print(3)

def func4(we):
    for oo in range(271,361):
        g(oo)
    print(4)

pool = multiprocessing.Pool(processes=4)        

pool.apply_async(func1, (1, ))
pool.apply_async(func2, (2, ))
pool.apply_async(func3, (3, ))
pool.apply_async(func4, (4, ))

pool.close()
pool.join()

end=time.time()
print(end-s1)

