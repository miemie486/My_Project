import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

f=open("data.txt")
t=f.readlines()
x=[]
y=[]
p=len(t)
print(p)
i=0
g=0
while(i<p-1):
    t1=t[i].split()
    x.append(float(t1[0]))
    y.append(float(t1[1]))
    g=g+3/800*float(t1[1])
    i=i+1
print(g)
plt.plot(x,y)
plt.show()