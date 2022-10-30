import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
import N_doping as N_D

def importfile(filename):
    f0=open(filename)
    f1=f0.readlines()
    x=[]
    y=[]
    for i in f1:
        t0=i.split()
        x.append(float(t0[0])*1e-6)
        y.append(float(t0[1])*1e6)
    return [x,y]

def Line_1(x):
    a=27.9175-3.35184*x
    return a
'''
def Line_2(x):
    a=25.4118+0.00858*x
    return a
'''
def Line_2(x):
    a=25.3358+0.06456*x-0.008811*x**2
    return a

data=importfile('data_800_20min.txt')

y_line_1=[]
y_line_2=[]
y_threshold=[]
iii=0
for i in data[0]:
    y_line_1.append(abs(Line_1(i*1e6)-m.log10(data[1][iii])))
    y_line_2.append(abs(Line_2(i*1e6)-m.log10(data[1][iii])))
    y_threshold.append(0.05)
    iii=iii+1

data[0]=np.array(data[0])*10**6

plt.figure(1)
plt.xlim(0,5)
plt.ylim(0,0.5)
plt.tick_params(labelsize=10)
#plt.semilogy(data[0],y_exp,label="Experiment Data")
plt.scatter(data[0],y_line_1,label="Reaction and Diffusion Points")
plt.scatter(data[0],y_line_2,label="Diffusion Points")
plt.plot(data[0],y_threshold,color='black')
plt.xlabel("depth, μm")
plt.ylabel("error, arb. unit")
plt.title("800℃ 20min 3.3Pa")
plt.legend(loc=0,ncol=1)
plt.show()