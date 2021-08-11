import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

f1=open("./论文氮扩散数据/[7]N_2min_3.3_5ep.txt")
f2=open("./论文氮扩散数据/[7]N_2min_6min_3.3_ep6.txt")
f3=open("N_2min_6min.txt")
f4=open("N_2min_6min_2.txt")

t1=f1.readlines()
t2=f2.readlines()
t3=f3.readlines()
t4=f4.readlines()
p=len(t1)

x1=[]
x2=[]
x3=[]
x4=[]
y1=[]
y2=[]
y3=[]
y4=[]

for i in range(p):
  x1.append(float(t1[i].split()[0])*1e6)
  x2.append(float(t2[i].split()[0])*1e6)
  x3.append(float(t3[i].split()[0]))
  x4.append(float(t4[i].split()[0]))

  y1.append((float(t1[i].split()[1])+10**23)/(5e28))
  y2.append((float(t2[i].split()[1])+10**23)/(5e28))
  y3.append(float(t3[i].split()[1])+10**23)
  y4.append(float(t4[i].split()[1])+10**23)
plt.figure(1)
plt.xlim(0,20)
plt.ylim(0,0.5)
plt.title("2N0 and 2N6 3.3Pa")
plt.xlabel("μm")
plt.ylabel("u/u_max")
plt.axvline(x=5,c="black",label="x=5μm")
#plt.axvline(x=6,c="red",label="x=6μm")
#plt.axvline(x=30,c="brown",label="x=30μm")
plt.plot(x1,y1,label="N distribution 2N0")
plt.plot(x2,y2,label="N distribution 2N6")
#plt.semilogy(x2,y2,label="one")
plt.legend(loc=0,ncol=1)

"""
plt.figure(2)
plt.xlim(0,30*10**-6)
plt.semilogy(x3,y3,label="two")
plt.semilogy(x4,y4,label="one")
plt.legend(loc=0,ncol=1)
"""

plt.show()
