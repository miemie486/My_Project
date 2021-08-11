import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

T=2
hbar=1.05457*10**(-34)
kb=1.38*10**(-23)
omega=1.3*10**9
miu0=4*m.pi*10**-7
Tc=9.23
delta_0=3.85*kb*Tc/2
rou_s=1.45*10**-7

def lambda_gamma(gamma):
  a=2*miu0*delta_0/(hbar*rou_s)*\
    (m.atan(delta_0/gamma)\
      -(m.pi**2*kb**2*T**2*gamma*delta_0)/(3*(gamma**2+delta_0**2)**2))
  return m.sqrt(1/a)

def lam(gamma,x):
  a=2*miu0*x*delta_0/(hbar*rou_s)*\
    (m.atan(x*delta_0/gamma)\
      -(m.pi**2*kb**2*T**2*gamma*x*delta_0)/(3*(gamma**2+(x*delta_0)**2)**2))
  #return 1/a
  return m.sqrt(1/a)
x=np.linspace(0.99,1,101)
y=np.linspace(0.00001,0.5,101)
X,Y=np.meshgrid(x,y)
Z=np.sin(X)+np.cos(Y)
z=[]
for i in x:
  y_r=[]
  for j in y:
    y_r.append(lam(i*delta_0,j))
  #print(y_r)
  z.append(y_r)

z=np.matrix(z)

fig = plt.figure()  #定义新的三维坐标轴
ax2 = Axes3D(fig)
ax3=plt.axes(projection='3d')
ax3.plot_surface(X,Y,z,cmap='rainbow')

plt.show()
