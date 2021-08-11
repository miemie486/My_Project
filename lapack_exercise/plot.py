import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

f=open("data.txt")

data0=f.readlines()

x=[]
y=[]
xx=np.linspace(0,5,1000)
yy=xx*0
for i in range(len(data0)):
    t0=data0[i].split()
    t1=list(map(float,t0))
    if(t1[1]>=40 or t1[1]<=-40):
        continue
    else:
        y.append(t1[1])
        x.append(t1[0])
plt.plot(x,y)
#plt.plot(xx,yy)
plt.show()