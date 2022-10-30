import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
#批量导入模块
k=['1','2','3','4']
i=0
f=[]
p=[]
t3=[]
tt3=[]
while i<len(k): 
    a='data'
    b='.txt'
    a=a+k[i]+b
    f.append(open(a))
    t=f[i].readline()
    t1=t.split()
    t2=list(map(float,t1))
    t2=np.array(t2)
    t3.append(t2)

    tt=f[i].readline()
    tt1=tt.split()
    tt2=list(map(float,tt1))
    p.append(len(tt2))
    tt2=np.array(tt2)
    tt3.append(tt2)
    i=i+1
print(t3)
print(tt3)
#总画图模块
n=2 #用于拟合的多项式阶数
bkb=0
while bkb<len(t3):
    i=0
    a1=np.zeros((n,n+1))
    while i<n:
        j=0
        while j<n+1:
            k=0
            while k<len(t3[bkb]):
                if j<n:
                    if i==0 and j==0:
                        a1[i][j]=n
                    else:
                        a1[i][j]=a1[i][j]+t3[bkb][k]**(i+j)
                if j==n:
                    a1[i][j]=a1[i][j]+tt3[bkb][k]*t3[bkb][k]**i
                k=k+1
            j=j+1
        i=i+1
    a1=np.array(a1)
    i=0
    while i<n:
        if i==0:
            a1[i]=a1[i]/a1[i][i]
        else:
            j=0
            while j<i:
                a1[i]=a1[i]-a1[j]*a1[i][j]/a1[j][j]
                a1[i]=a1[i]/a1[i][i]
                j=j+1
        i=i+1
    i=n-1
    while i>=0:
        j=n-1
        while j>i:
            a1[i]=a1[i]-a1[j]*a1[i][j]
            j=j-1
        i=i-1
    plt.figure(bkb)
    x=np.linspace(5,13,1000)
    y=a1[0][n]+a1[1][n]*x
    plt.plot(x,y)
    plt.scatter(t3[bkb],tt3[bkb])  # Only the first fit is good fit. Other fits are bad fit.
    bkb=bkb+1
plt.show()
