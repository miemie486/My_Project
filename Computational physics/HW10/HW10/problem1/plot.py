import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
f=open('data.txt')
t=f.readline()
t1=t.split()
t2=list(map(float,t1))
t2=np.array(t2)

tt=f.readline()
tt1=tt.split()
tt2=list(map(float,tt1))
tt2=np.array(tt2)

plt.plot(t2,tt2)
plt.show()