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

print(p)
h=[]
i=0
while i<=p-2:
    a=t3[i+1]-t3[i]
    h.append(a)
    i=i+1
i=0
#构建增广矩阵
matrix=[]
while i<=p-1:
    if i==0:
        j=1
        row=[1.0]
        while j<=p-1:
            row.append(0.0)
            j=j+1
        row.append(0.0)
        matrix.append(row)
        
    if i==1:
        j=3
        row=[0.0,2*(h[1]+h[2]),h[2]]
        while j<=p-1:
            row.append(0.0)
            j=j+1
        a=(1/h[i])*(tt3[i+1]-tt3[i])-(1/h[i-1])*(tt3[i]-tt3[i-1])
        a=6*a
        row.append(a)
        matrix.append(row)

    if (i>1)&(i<p-2):
        j=0
        row=[]
        while j<=p-1:
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
        a=(1/h[i])*(tt3[i+1]-tt3[i])-(1/h[i-1])*(tt3[i]-tt3[i-1])
        a=6*a
        row.append(a)
        matrix.append(row)
    if i==p-2:
        j=0
        row=[]
        while j<=p-1:
            if j<=p-4:
                row.append(0.0)
                j=j+1
            else:
                row.append(h[i-1])
                row.append(2*(h[i-1]+h[i]))
                row.append(0)
                j=j+3
        a=(1/h[i])*(tt3[i+1]-tt3[i])-(1/h[i-1])*(tt3[i]-tt3[i-1])
        a=6*a
        row.append(a)
        matrix.append(row)
    if i==p-1:
        j=0
        row=[]
        while j<=p-1:
            if j!=p-1:
                row.append(0.0)
                j=j+1
            else:
                row.append(1.0)
                j=j+1
        a=0.0
        row.append(a)
        matrix.append(row)
    i=i+1

Matrix=np.array(matrix)

#输出未化为上对角阵的原矩阵
kk=0
while kk<=p-1:
    bb=0
    while bb<=p:
        a=Matrix[kk][bb]
        if bb!=p:
            print("%.3f"%a,end=" ")
        else:
            print("%.3f"%a)
        bb=bb+1
    kk=kk+1
print("______________")
#高斯消元法
i=1
while i<=p-2:
    if i<=p-3:
        Matrix[i]=Matrix[i]/Matrix[i][i]
        Matrix[i+1]=Matrix[i+1]-Matrix[i]*Matrix[i+1][i]
        i=i+1
    else:
        Matrix[i]=Matrix[i]/Matrix[i][i]
        i=i+1
i=p-2
while i>=1:
    Matrix[i-1]=Matrix[i-1]-Matrix[i]*Matrix[i-1][i]
    i=i-1

#输出化为对角阵后的矩阵
kk=0
while kk<=p-1:
    bb=0
    while bb<=p:
        a=Matrix[kk][bb]
        if bb!=p:
            print("%.3f"%a,end=" ")
        else:
            print("%.3f"%a)
        bb=bb+1
    kk=kk+1

def pp(j):
    a=Matrix[j][p]
    return a

def p_(x,j):
    j=j-1
    a=((pp(j+1)-pp(j))/(6*h[j]))*(x-t3[j])**3 \
        +(pp(j)/2)*(x-t3[j])**2 \
        +((tt3[j+1]-tt3[j])/h[j]-(h[j]*pp(j+1))/6-(h[j]*pp(j))/3)*(x-t3[j]) \
        +tt3[j]
    return a
#检验插值函数的一阶导是否连续
j=0
while j<=p-2:
    kkk=(tt3[j+1]-tt3[j])/h[j]-(h[j]*pp(j+1))/6-h[j]*pp(j)/3
    kkkk=kkk+(pp(j+1)-pp(j))*h[j]**2/(2*h[j])+pp(j)*h[j]
    print("p'(",j,")=",kkk)
    print("p'(",j+1,")=",kkkk)
    j=j+1
    
j=1
xx=[]
yy=[]
while j<=p-1:
    k=1001     #每两个点之间的插值点数
    x=np.linspace(t3[j-1],t3[j],k)
    i=0
    y=[]
    while i<=k-1:
        y.append(p_(x[i],j))
        i=i+1
    yy=yy+y
    xx=np.append(xx,x)
    j=j+1

plt.plot(xx,yy)
plt.scatter(t3,tt3)
plt.show()

    


