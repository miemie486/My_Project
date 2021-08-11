import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

T=2
hbar=1.05457*10**(-34)
kb=1.38*10**(-23)
omega=2*m.pi*1.3*10**9
lambda0=32*10**-9
xi_0=40*10**-9
l1=3*10**-9
l2=20*10**-9
ds=50*10**-9


def lambda_l(l):
  a=lambda0*(1+xi_0/l)**0.5
  return a

def B1(x):
  a=(m.cosh((ds-x)/lambda_l(l1))+(lambda_l(l2)/lambda_l(l1))*m.sinh((ds-x)/lambda_l(l1)))\
    /(m.cosh(ds/lambda_l(l1))+lambda_l(l2)/lambda_l(l1)*m.sinh(ds/lambda_l(l1)))
  return a
def B3(x):
  a=m.exp(-(x-ds)/lambda_l(l2))/(m.cosh(ds/lambda_l(l1))+lambda_l(l2)/lambda_l(l1)*m.sinh(ds/lambda_l(l1)))
  return a
x=np.linspace(0,ds,100001)
y=[]
y1=[]
yy=[]
yy1=[]
for i in x:
  y.append(B1(i))
  y1.append(m.exp(-i/lambda_l(l1)))
xx=np.linspace(ds,3*ds,100001)
for i in xx:
  yy.append(B3(i))
  yy1.append(m.exp(-(i-ds)/lambda_l(l2))*m.exp(-ds/lambda_l(l1)))
#print((f(0.5*xi+0.0001*xi)-f(0.5*xi-0.0001*xi))/(0.0002*xi))
plt.plot(x,y)
plt.plot(xx,yy)
plt.plot(x,y1)
plt.plot(xx,yy1)
#plt.plot(yy,y)
plt.show()
