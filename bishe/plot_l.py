import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

T=2
hbar=1.05457*10**(-34)
kb=1.38*10**(-23)
omega=2*m.pi*1.3*10**9




def f(l):
  aa=2*l/(m.sqrt(2)*kb*T)*(hbar*omega/(m.pi*l))**1.5*m.log(l/(hbar*omega))*m.exp(-l/(kb*T))/(m.tanh(l/(2*kb*T))**1.5*(1+m.exp(-l/kb*T))**2)
  return aa


x=np.linspace(0.999,0.02,100001)
y=[]
yy=[]
for i in x:
  print(i)
  y.append(f(i*1.765*kb*9.23))
  #yy.append(0.5*xi)
#print((f(0.5*xi+0.0001*xi)-f(0.5*xi-0.0001*xi))/(0.0002*xi))
plt.semilogy(x,y)
#plt.plot(yy,y)
plt.show()
