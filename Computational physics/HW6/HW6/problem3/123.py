import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

list1=[5,7,9,17,19,21]

def f(x):
    a=x*m.sin(2*m.pi*x+1)
    return a

def L(i,n,x,nn):#nn表示第nn个的系数
    global x2
    k=1
    j=0
    while j<=n-1:
        if i!=j:
            k=k*((x-x2[nn][j])/(x2[nn][i]-x2[nn][j]))
        j=j+1
    return k
def LL(x,nn):
    kk=0
    a=0
    while kk<=list1[nn]-1:
        a=a+L(kk,list1[nn],x,nn)*y2[nn][kk]
        kk=kk+1
    return a


x1=np.linspace(-1,1,1001)
p1=len(x1)
y1=[]
i=0
while i<=p1-1:
    y1.append(f(x1[i]))
    i=i+1



i=0
x2=[]
y2=[]
while i<=len(list1)-1:

    x2.append(np.linspace(-1,1,list1[i]))
    p2=len(x2[i])
    
    y22=[]
    j=0
    while j<=p2-1:
        y22.append(f(x2[i][j]))
        j=j+1
    y2.append(y22)
    i=i+1

iii=0
while iii<=len(list1)-1:
    xx=np.linspace(-1,1,1001)
    yy=[]
    ii=0
    while ii<=1000:
        yy.append(LL(xx[ii],iii))
        ii=ii+1
    yy=np.array(yy)
    plt.plot(xx,yy)
    plt.scatter(x2[iii],y2[iii])
    iii=iii+1

plt.plot(x1,y1)
plt.ylim(-1.5,1.5)
plt.show()