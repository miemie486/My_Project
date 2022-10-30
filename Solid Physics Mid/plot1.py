import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
a=1
def f(x,y):
    oo=3.0-m.cos((1.0/m.sqrt(3.0))*y*a-x*a)-m.cos(a*2.0*y/m.sqrt(3.0))-m.cos((-1.0/m.sqrt(3.0))*y*a-x*a)
    return 2*oo

x=np.linspace(-5*m.pi,5*m.pi,500)
y=[]
for i in x:
    y.append(m.sqrt(f(i,0)))
plt.plot(x,y)
plt.show()