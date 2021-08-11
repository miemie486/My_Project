import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

T=2
hbar=1.05457*10**(-34)
kb=1.38*10**(-23)
omega=1.3*10**9

t=3600*3

D0=2.6*10**-6

D=D0*m.exp(-(1.523*10**5/(6.02*10**23))/(kb*(800+273)))
D=1.37*10**-13
print(D)
def f(l):
  aa=m.erfc(l/(2*(D*t)**0.5))
  #print((2*(D*t)**0.5))
  return aa


x=np.linspace(0,400*10**-6,10000)
y=[]
yy=[]
for i in x:
  y.append(f(i))
  #yy.append(0.5*xi)
#print((f(0.5*xi+0.0001*xi)-f(0.5*xi-0.0001*xi))/(0.0002*xi))
plt.semilogy(x,y)
#aaa=[900,950,1000,1100,1200,1300]
#bbb=[0.8*10**-7,1.53*10**-7,2.25*10**-7,4.14*10**-7,9.86*10**-7,1.2*10**-6]
aaa=[900,950,1000,1100,1200]
bbb=[0.8*10**-7,1.53*10**-7,2.25*10**-7,4.14*10**-7,9.86*10**-7]

"""

omg=np.linspace(22.9,23,101)
for omg_1 in omg:
  ccc=[]
  fangcha=0
  j=0
  for i in aaa:
    ccc.append(D0*m.exp(-(22.9*4.186*10**3/(6.02*10**23))/(kb*(i+273))))
    fangcha=fangcha+(bbb[j]-ccc[j])**2
    j=j+1
  biaozhuncha=(fangcha/6)**0.5
  #print(biaozhuncha*10**7,omg_1)
plt.scatter(aaa,bbb)
plt.plot(aaa,ccc)

"""

#plt.plot(yy,y)
plt.show()
