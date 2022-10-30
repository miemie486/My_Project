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

#f1=np.load("data_without_transition_region\data_800_20min_theory_with_transition_regin.npy")
#f1=np.load("data\_0min baking\data_800_3N60min_theory.npy")
#f1=np.load("data\data_800_3N60min_theory.npy")
f1=np.load("data\_0min baking\data_800_10N0min_theory.npy")
#f1=np.load("data\\variable_pressure_N2A6\\data_800_2N6_min_10.0Pa_theory.npy")
#f1=np.load("1_test_data_800_20N0min_theory.npy")
#f1=np.load("_test_data_800_3N60min_theory.npy")
data=importfile('data_800_20min.txt')

X=np.linspace(0,N_D.thickness_a,400)
print(np.interp(15*10**-6,X,f1))

X=X*10**6
data[0]=np.array(data[0])*10**6

plt.figure(1)
plt.xlim(0,15)
plt.ylim(1e24,1e28)
#plt.ylim(1e25,1e26)
plt.tick_params(labelsize=10)
plt.semilogy(X,f1+1e18,label="Theory_1")
#plt.semilogy(X,f2+1e18,label="Theory_2")
plt.semilogy(data[0],data[1],label="Experiment")
plt.xlabel("depth, μm")
plt.ylabel("n, atom/m$^3$")
plt.title("800℃ 20min 3.3Pa")
plt.legend(loc=0,ncol=1)
#plt.show()

