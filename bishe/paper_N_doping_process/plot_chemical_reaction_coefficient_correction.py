
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
    a=25.3358+0.06456*x-0.0088109*x**2
    return a

data=importfile('data_800_20min.txt')

y_line_1=[]
y_line_2=[]
y_exp=[]

for i in data[0]:
    y_line_1.append(10**Line_1(i*1e6))
    y_line_2.append(10**Line_2(i*1e6))

data[0]=np.array(data[0])*10**6

divide1_2=np.ones(len(data[0]))*0.600
divide2_3=np.ones(len(data[0]))*0.934

font={}


plt.figure(1)
plt.xlim(0,5)
plt.ylim(1e24,1e28)
plt.tick_params(labelsize=10)
plt.semilogy(data[0],data[1],label="Experiment Data")
#plt.semilogy(data[0],y_line_1,label="Reaction and Diffusion Line")
#plt.semilogy(data[0],y_line_2,label="Diffusion Line")
#font=
plt.text(0.18,5*1e27,'(1)',fontsize=12,color='r')
plt.text(0.72,5*1e27,'(2)',fontsize=12,color='r')
plt.text(1.5,5*1e27,'(3)',fontsize=12,color='r')

plt.vlines(0.56,1e28,1e21,colors="black")
plt.vlines(1.10,1e28,1e21,colors="black")
plt.xlabel("depth, μm")
plt.ylabel("n, atom/m$^3$")
plt.title("800℃ 20min 3.3Pa")
plt.legend(loc=0,ncol=1)
plt.show()
