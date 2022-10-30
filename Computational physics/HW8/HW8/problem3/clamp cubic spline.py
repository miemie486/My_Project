import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

def fg(x):
    aa=x*m.sin(2*m.pi*x+1)
    return aa

def Df(x):
    hh=0.00000001
    aaa=(fg(x-2*hh)-8*fg(x-hh)+8*fg(x+hh)-fg(x+2*hh))/(12*hh)
    return aaa

#批量导入模块
k=['9','17','19','21']
lenk=len(k)
kk123=['9','17','19','21','y=x*sin(2*pi*x+1)','9','17','19','21']
i=0
f=[]
p=[]
t3=[]
tt3=[]
while i<lenk: 
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

#总画图模块
bkb=0
while bkb<lenk:
    h=[]
    i=0
    while i<=p[bkb]-2:
        a=t3[bkb][i+1]-t3[bkb][i]
        h.append(a)
        i=i+1
    i=0
    matrix=[]
    while i<=p[bkb]-1:
        if i==0:
            j=2
            row=[2*h[0],h[0]]
            while j<=p[bkb]-1:
                row.append(0.0)
                j=j+1
            afk=t3[bkb][0]
            print(afk)
            a=(1/h[i])*(tt3[bkb][i+1]-tt3[bkb][i])-Df(afk) 
            a=6*a
            row.append(a)
            matrix.append(row)

        if i==1:
            j=3
            row=[h[0],2*(h[0]+h[1]),h[1]]
            while j<=p[bkb]-1:
                row.append(0.0)
                j=j+1
            a=(1/h[i])*(tt3[bkb][i+1]-tt3[bkb][i])-(1/h[i-1])*(tt3[bkb][i]-tt3[bkb][i-1])
            a=6*a
            row.append(a)
            matrix.append(row)

        if (i>1)&(i<p[bkb]-1):
            j=0
            row=[]
            while j<=p[bkb]-1:
                if (j!=i)&(j!=i-1)&(j!=i+1):
                    row.append(0.0)
                    j=j+1
                else:
                    j=j+1
                if j==i-1:
                    row.append(h[i-2])
                    row.append(2*(h[i-2]+h[i-1]))
                    row.append(h[i-1])
                    j=j+3
            a=(1/h[i])*(tt3[bkb][i+1]-tt3[bkb][i])-(1/h[i-1])*(tt3[bkb][i]-tt3[bkb][i-1])
            a=6*a
            row.append(a)
            matrix.append(row)

        if i==p[bkb]-1:
            j=0
            row=[]
            while j<=p[bkb]-2:
                if j!=p[bkb]-2:
                    row.append(0.0)
                    j=j+1
                else:
                    row.append(h[j])
                    row.append(2*h[j])
                    j=j+1
            a=Df(t3[bkb][i])-(1/h[i-1])*(tt3[bkb][i]-tt3[bkb][i-1])
            a=6*a
            row.append(a)
            matrix.append(row)
        i=i+1
    print(matrix)
    Matrix=np.array(matrix)

    #输出未化为上对角阵的原矩阵
    kk=0
    while kk<=p[bkb]-1:
        bb=0
        while bb<=p[bkb]:
            a=Matrix[kk][bb]
            if bb!=p[bkb]:
                print("%.3f"%a,end=" ")
            else:
                print("%.3f"%a)
            bb=bb+1
        kk=kk+1
    print("______________")
    #高斯消元法
    i=0
    while i<=p[bkb]-1:
        if i<=p[bkb]-2:
            Matrix[i]=Matrix[i]/Matrix[i][i]
            Matrix[i+1]=Matrix[i+1]-Matrix[i]*Matrix[i+1][i]
            i=i+1
        else:
            Matrix[i]=Matrix[i]/Matrix[i][i]
            i=i+1
    i=p[bkb]-1
    while i>=0:
        Matrix[i-1]=Matrix[i-1]-Matrix[i]*Matrix[i-1][i]
        i=i-1

    #输出化为对角阵后的矩阵
    kk=0
    while kk<=p[bkb]-1:
        bb=0
        while bb<=p[bkb]:
            a=Matrix[kk][bb]
            if bb!=p[bkb]:
                print("%.3f"%a,end=" ")
            else:
                print("%.3f"%a)
            bb=bb+1
        kk=kk+1

    def pp(j):
        a=Matrix[j][p[bkb]]
        return a

    def p_(x,j):
        j=j-1
        a=((pp(j+1)-pp(j))/(6*h[j]))*(x-t3[bkb][j])**3 \
            +(pp(j)/2)*(x-t3[bkb][j])**2 \
            +((tt3[bkb][j+1]-tt3[bkb][j])/h[j]-(h[j]*pp(j+1))/6-(h[j]*pp(j))/3)*(x-t3[bkb][j]) \
            +tt3[bkb][j]
        return a
    #检测插值函数的一阶导是否连续
    j=0
    while j<=p[bkb]-2:
        kkk=(tt3[bkb][j+1]-tt3[bkb][j])/h[j]-(h[j]*pp(j+1))/6-h[j]*pp(j)/3
        kkkk=kkk+(pp(j+1)-pp(j))*h[j]**2/(2*h[j])+pp(j)*h[j]
        print("p'(",j,")=",kkk)
        print("p'(",j+1,")=",kkkk)
        j=j+1
    #构建画图所需数组并画图
    j=1
    xx=[]
    yy=[]
    while j<=p[bkb]-1:
        k=101
        x=np.linspace(t3[bkb][j-1],t3[bkb][j],k)
        i=0
        y=[]
        while i<=k-1:
            y.append(p_(x[i],j))
            i=i+1
        yy=yy+y
        xx=np.append(xx,x)
        j=j+1
    plt.plot(xx,yy)
    plt.scatter(t3[bkb],tt3[bkb])
    bkb=bkb+1
xxx=np.linspace(-1,1,1000)
yyy=xxx*np.sin(2*np.pi*xxx+1)
plt.plot(xxx,yyy)
plt.legend(kk123)
plt.show()
