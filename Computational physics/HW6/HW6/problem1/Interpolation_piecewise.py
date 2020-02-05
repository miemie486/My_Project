import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
f=open("777.txt")

t=f.readline()
t1=t.split()
t2=list(map(float,t1))
t3=np.array(t2)

tt=f.readline()
tt1=tt.split()
tt2=list(map(float,tt1))
tt3=np.array(tt2)
p=len(tt3)

def f(i,j,n,x):#从第j个点往后数n个点的拉格朗日插值多项式的从j往后数第i个点的系数
    global t3
    k=1
    jj=j-1
    j=j-1
    i=i+jj
    l=0
    while j<=n-1+jj:
        if i!=j:
            k=k*((x-t3[j])/(t3[i]-t3[j]))
            l=l+1
        j=j+1
    return k
def ff(x,j,n):#从第j个点往后数n个点所拟合的n-1次拉格朗日插值多项式
    kk=0
    a=0
    while kk<=n-1:
        a=a+f(kk,j,n,x)*tt3[kk+j-1]
        kk=kk+1
    return a
yy=[]
xxx=[]
n=5 #n-1为分段插值所用的插值多项式阶数
gg=0
nn=1001 #每n个数据点之间画的点的个数
while gg<=p-n:
    if gg==0:
        xx=[]
        xx=np.linspace(t3[gg],t3[gg+1],nn)
        xxx=np.append(xxx,xx)
        ii=0
        while ii<=nn-1:
            yy.append(ff(xx[ii],gg+1,n))
            
            ii=ii+1
        gg=gg+1

    xx=[]
    xx=np.linspace(t3[gg],t3[gg+n-1],nn)
    xxx=np.append(xxx,xx)
    ii=0
    while ii<=nn-1:
        yy.append(ff(xx[ii],gg+1,n))  
        ii=ii+1
    gg=gg+1

    if gg==p-n:
        xx=[]
        xx=np.linspace(t3[gg+n-2],t3[gg+n-1],nn)
        xxx=np.append(xxx,xx)
        ii=0
        while ii<=nn-1:
            yy.append(ff(xx[ii],gg+1,n))
            
            ii=ii+1
        gg=gg+1

yy=np.array(yy)
plt.plot(xxx,yy)
plt.scatter(t3,tt3)
plt.show()



