import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

f1=open("N_2min_3.3.txt")
f2=open("N_2min_2.txt")
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
  x1.append(float(t1[i].split()[0]))
  x2.append(float(t2[i].split()[0]))
  x3.append(float(t3[i].split()[0]))
  x4.append(float(t4[i].split()[0]))

  y1.append(float(t1[i].split()[1])+10**23)
  y2.append(float(t2[i].split()[1])+10**23)
  y3.append(float(t3[i].split()[1])+10**23)
  y4.append(float(t4[i].split()[1])+10**23)
plt.figure(1)
plt.xlim(0,10*10**-6)
plt.semilogy(x1,y1,label="1")
#plt.semilogy(x2,y2,label="one")
plt.legend(loc=0,ncol=1)

plt.figure(2)
plt.xlim(0,30*10**-6)
plt.semilogy(x3,y3,label="two")
#plt.semilogy(x4,y4,label="one")
plt.legend(loc=0,ncol=1)

plt.show()
