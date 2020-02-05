import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
import multiprocessing
import csv
import time
from numba import jit
#内存为16G的电脑用本程序画出来的分型图的极限是15001*15001

s1=time.time()
@jit
def fractal(xx,ymin,ymax,precision):#xx为x轴上数据点个数，从y轴上第ymin个数据点开始，第ymax个数据点结束，precision是精度
    global x,y
    z=[]
    ky=ymin
    while ky<ymax and ky>=ymin:
        kx=0
        zx=[]
        while kx<=xx:
            i=0
            zi=x[kx]+y[ky]*1j
            zi1=zi
            zi_plus=zi+1
            while abs(zi1-zi_plus)>=precision and zi!=0:
                zi_plus=zi-(zi**3-1)/(3*zi**2)
                zi1=zi
                zi=zi_plus
                i=i+1
            zx.append(i)
            kx=kx+1
        ky=ky+1
        z.append(zx)
    return z

jingdu=0.0001   #精度
x=np.linspace(-2.5,2.5,15001)
y=np.linspace(-2.5,2.5,15001)
px=len(x)-1
p=int(len(y)/4)
X,Y=np.meshgrid(x,y)

def func1(we):
    z1=fractal(px,0,p,jingdu)
    z1=np.array(z1)
    np.save("z1.npy",z1)
    print("233333")

def func2(we):
    z2=fractal(px,p,p+p,jingdu)
    z2=np.array(z2)
    np.save("z2.npy",z2)
    print("233333")

def func3(we):
    z3=fractal(px,p+p,p+2*p,jingdu)
    z3=np.array(z3)
    np.save("z3.npy",z3)
    print("233333")

def func4(we):
    z4=fractal(px,p+2*p,p+3*p+1,jingdu)
    z4=np.array(z4)
    np.save("z4.npy",z4)
    print("233333")
if __name__ == "__main__":
    pool = multiprocessing.Pool(processes=4)        

    pool.apply_async(func1, (1, ))
    pool.apply_async(func2, (2, ))
    pool.apply_async(func3, (3, ))
    pool.apply_async(func4, (4, ))

    pool.close()
    pool.join()


s2=time.time()


z1=np.load("z1.npy")
z1=list(z1)

z2=np.load("z2.npy")
z2=list(z2)

z3=np.load("z3.npy")
z3=list(z3)

z4=np.load("z4.npy")
z4=list(z4)

s3=time.time()

z=z1+z2+z3+z4
z=np.array(z)

plt.figure(figsize=(150,150))
plt.contourf(X,Y,z)
plotPath='fractal996.png'
#plt.savefig(plotPath)
#plt.show()

end=time.time()
p1=(s2-s1)
p2=(s3-s2)
p3=(end-s3)
p4=(end-s1)
print("算点耗时：",p1)
print("IO耗时：",p2)
print("画图耗时：",p3)
print("总耗时",p4)