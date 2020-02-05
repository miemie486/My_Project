import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.figure(1)
name=['1','2','3']
label_=['Euler','Midpoint','Heuns']
xx=[]
yy=[]
for i in range(0,3):
    a1=name[i]
    a2='.txt'
    a3=a1+a2
    b1=open(a3)
    b2=b1.read()
    c1=b2.split()
    c2=list(map(float,c1))
    d1=[]
    d2=[]
    j=0
    while j<len(c2):
        d1.append(10**-c2[j])
        d2.append(c2[j+1])
        j=j+2
    print(d1)
    print(d2)
    plt.loglog(d1,d2,label=label_[i])
    plt.xlabel('h')
    plt.ylabel('error')
    plt.legend()
plt.show()
#如果误差下降到10^-12以下之后，就需要考虑舍入误差。
    