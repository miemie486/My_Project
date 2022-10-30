from cmath import inf
import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
import N_doping as N_D

#判定标准见文章
def f(x,X,f1):
    a=np.interp(x,X,f1)-1.09e25#1.14e25
    return a

def f_top(x,X,f1):
    a=np.interp(x,X,f1)-1.62e25#1.85e25
    return a

def f_bottom(x,X,f1):
    a=np.interp(x,X,f1)-7.25e24#7.07e24
    return a    
       
    
baking_time="6"
doping_time=[2]
doping_time.sort()

pressure_list=[]
for i in range(1,101):
    pressure_list.append(i/10)
print(pressure_list)


polish_depth_best=[]
for j in pressure_list:
    for i in doping_time:
        name="data\\variable_pressure_N"+str(i)+"A"+baking_time+"\\data_800_"+str(i)+"N"+baking_time+"_min_"+str(j)+"Pa_theory.npy"
        f1=np.load(name)
        X=np.load("data\X.npy")
        
        xn=5*10**-6
        xn_plus_1=6*10**-6
        error=abs(xn-xn_plus_1)*10**6
        while(error>=0.01):
            xn_plus_1_local=xn_plus_1
            xn_plus_1=xn_plus_1-f(xn_plus_1,X,f1)*(xn_plus_1-xn)/(f(xn_plus_1,X,f1)-f(xn,X,f1))
            xn=xn_plus_1_local
            error=abs(f(xn_plus_1,X,f1)*1e-25)

        if(m.isnan(xn_plus_1)):
            xn_plus_1=0
        polish_depth_best.append(xn_plus_1)
        print(i,xn_plus_1,m.isnan(xn_plus_1))


polish_depth_top=[]
for j in pressure_list:
    for i in doping_time:
        name="data\\variable_pressure_N"+str(i)+"A"+baking_time+"\\data_800_"+str(i)+"N"+baking_time+"_min_"+str(j)+"Pa_theory.npy"
        f1=np.load(name)
        X=np.load("data\X.npy")
        
        
        xn=5*10**-6
        xn_plus_1=6*10**-6
        error=abs(xn-xn_plus_1)*10**6

        #print(f_top(xn_plus_1,X,f1)-f_top(xn,X,f1))

        while(error>=0.01):
            xn_plus_1_local=xn_plus_1
            xn_plus_1=xn_plus_1-f_top(xn_plus_1,X,f1)*(xn_plus_1-xn)/(f_top(xn_plus_1,X,f1)-f_top(xn,X,f1))
            xn=xn_plus_1_local
            error=abs(f_top(xn_plus_1,X,f1)*1e-25)

        if(m.isnan(xn_plus_1)):
            xn_plus_1=0
        polish_depth_top.append(xn_plus_1)
        print(i,xn_plus_1,m.isnan(xn_plus_1))

polish_depth_bottom=[]
for j in pressure_list:
    for i in doping_time:
        name="data\\variable_pressure_N"+str(i)+"A"+baking_time+"\\data_800_"+str(i)+"N"+baking_time+"_min_"+str(j)+"Pa_theory.npy"
        f1=np.load(name)
        X=np.load("data\X.npy")
        
        xn=5*10**-6
        xn_plus_1=6*10**-6
        error=abs(xn-xn_plus_1)*10**6
        while(error>=0.01):
            xn_plus_1_local=xn_plus_1
            xn_plus_1=xn_plus_1-f_bottom(xn_plus_1,X,f1)*(xn_plus_1-xn)/(f_bottom(xn_plus_1,X,f1)-f_bottom(xn,X,f1))
            xn=xn_plus_1_local
            error=abs(f_bottom(xn_plus_1,X,f1)*1e-25)

        if(m.isnan(xn_plus_1)):
            xn_plus_1=0
        polish_depth_bottom.append(xn_plus_1)
        print(i,xn_plus_1,m.isnan(xn_plus_1))


polish_depth_best=np.array(polish_depth_best)*10**6
polish_depth_top=np.array(polish_depth_top)*10**6
polish_depth_bottom=np.array(polish_depth_bottom)*10**6

plt.figure(1)
plt.tick_params(labelsize=10)
plt.plot(pressure_list,polish_depth_best,label="mathematical expectation")
plt.plot(pressure_list,polish_depth_top,label="upper limit of the interval")
plt.plot(pressure_list,polish_depth_bottom,label="lower limit of the interval")
plt.fill_between(pressure_list, polish_depth_top,polish_depth_bottom,facecolor='green', alpha=0.3)
plt.xlabel("nitrogen pressure, Pa")
plt.ylabel("polishing depth, μm")
plt.title("800℃ 3.3Pa "+baking_time+" annealing")
plt.legend(loc=0,ncol=1)
plt.show()

'''
plt.figure(1)
plt.xlim(0,5*10**-6)
plt.ylim(1e24,1e28)
plt.tick_params(labelsize=10)
plt.semilogy(X,f1+1e18,label="Theory")
plt.xlabel("depth, m")
plt.ylabel("n, atom/m^3")
plt.title("800_20min 3.3Pa")
plt.legend(loc=0,ncol=1)
plt.show()
'''