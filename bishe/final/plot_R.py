import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

f0=open("RBCS.txt")
t0=f0.readlines()
p=len(t0)    

t02x=[]
t02y=[]
for i in range(p-1):
  t02x.append(float(t0[i].split()[0])*1e9)
  t02y.append(float(t0[i].split()[1]))


plt.figure(1)
plt.loglog(t02x,t02y)
plt.grid()
plt.xlabel("Mean Free Path l [nm]")
plt.ylabel("R_BCS [Ω]")
plt.title("The Relationship Between R_BCS and Mean Free Path ")

plt.show()
