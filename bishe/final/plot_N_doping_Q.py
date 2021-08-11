import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

f0=open("droping.txt")
f1=open("data1.txt")
f2=open("./data/R_BCS_vol3/good_fit_data1.txt")

t0=f0.readlines()
t1=f1.readlines()
t2=f2.readlines()
p=len(t0)                   
t02x=[]
t02y=[]
t12x=[]
t12y=[]
t22x=[]
t22y=[]
for i in range(p-1):
  t02x.append(float(t0[i].split()[0])*4.25)
  t02y.append(float(t0[i].split()[1]))
  t12x.append(float(t1[i].split()[0]))
  t12y.append(float(t1[i].split()[1]))
  t22x.append(float(t2[i].split()[0]))
  t22y.append(float(t2[i].split()[1]))
plt.figure(1)

plt.axes(yscale = "log") 
plt.title("fit result")
plt.xlabel("B [mT]")
plt.ylabel("Q")
plt.scatter(t02x,t02y,c="black",label="experimental data")
#plt.plot(t12x,t12y)
plt.ylim(10**10,10**11)
plt.plot(t22x,t22y,label="fit curve")
plt.show()