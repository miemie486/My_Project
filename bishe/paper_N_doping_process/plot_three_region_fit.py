
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

font={}

region1_x=data[0][0:19]
region1_y=data[1][0:19]


region2_x=data[0][34:71]
region2_y=data[1][34:71]


plt.figure(1)
plt.xlim(0,0.5)
plt.ylim(1e26,1e28)
plt.tick_params(labelsize=10)
plt.semilogy(data[0],y_line_1,label="Reaction and Diffusion Line")
plt.scatter(region1_x,region1_y,label="Experiment Data",color="r")
plt.text(0.08,3*1e26,'y=10^(27.9175-3.3518x)',fontsize=13,color='black')
plt.xlabel("depth, μm")
plt.ylabel("n, atom/m$^3$")
plt.title("800℃ 20min 3.3Pa")
plt.legend(loc=0,ncol=1)

plt.figure(2)
plt.xlim(0.9,5.5)
plt.ylim(2.5e25,3e25)
plt.tick_params(labelsize=10)
plt.scatter(region2_x,region2_y,label="Experiment Data",color="r")
plt.text(2.3,2.6e25,'y=10^(25.3358+0.0646x-0.009x$^2$)',fontsize=13,color='black')
plt.semilogy(data[0],y_line_2,label="Diffusion Line")
plt.ylabel("n, atom/m$^3$")
plt.xlabel("depth, μm")
plt.title("800℃ 20min 3.3Pa")
plt.legend(loc=0,ncol=1)
plt.show()
