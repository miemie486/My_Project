from cmath import inf
import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
import N_doping as N_D
import os
from mpl_toolkits.mplot3d import Axes3D


#以四个工艺的最佳抛光深度对应的浓度的平均作为最佳浓度
#  2N0    5μm EP 0.98e25  
#  2N6    5μm EP 1.05e25 
#  20N0   12μm EP 1.14e25  
#  20N30  20μm EP 0.78e25
#  (0.983+1.05+1.14+0.782)/4=0.99e25  上下浮动20%作为冗余区间
def f(x,X,f1):
    a=np.interp(x,X,f1)-1e25
    return a

def f_top(x,X,f1):
    a=np.interp(x,X,f1)-1.2e25
    return a

def f_bottom(x,X,f1):
    a=np.interp(x,X,f1)-8e24
    return a    
    
name_list=os.listdir("data_3D")
doping_time=[]
baking_time=[]
doping_time_sort=[]
baking_time_sort=[]
for i in name_list:
    a=i.replace('data_800_','')
    a=a.replace('min_theory.npy','')
    b=a.split('N')
    doping_time.append(float(b[0]))
    baking_time.append(float(b[1]))
    doping_time_sort.append(float(b[0]))
    baking_time_sort.append(float(b[1]))


polish_depth_best=[]
for i in name_list:
    name="data_3D\\"+i

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
    #print(i,xn_plus_1,m.isnan(xn_plus_1))

polish_depth_top=[]
for i in name_list:
    name="data_3D\\"+i

    f1=np.load(name)
    X=np.load("data\X.npy")
    
    xn=5*10**-6
    xn_plus_1=6*10**-6
    error=abs(xn-xn_plus_1)*10**6
    while(error>=0.01):
        xn_plus_1_local=xn_plus_1
        xn_plus_1=xn_plus_1-f_top(xn_plus_1,X,f1)*(xn_plus_1-xn)/(f_top(xn_plus_1,X,f1)-f_top(xn,X,f1))
        xn=xn_plus_1_local
        error=abs(f_top(xn_plus_1,X,f1)*1e-25)

    if(m.isnan(xn_plus_1)):
        xn_plus_1=0
    polish_depth_top.append(xn_plus_1)
    #print(i,xn_plus_1,m.isnan(xn_plus_1))

polish_depth_bottom=[]
for i in name_list:
    name="data_3D\\"+i
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
    #print(i,xn_plus_1,m.isnan(xn_plus_1))

doping_time_sort.sort()
baking_time_sort.sort()

polish_depth_best=np.array(polish_depth_best)*10**6
polish_depth_top=np.array(polish_depth_top)*10**6
polish_depth_bottom=np.array(polish_depth_bottom)*10**6
doping_time=np.array(doping_time)
baking_time=np.array(baking_time)

Z_best=np.zeros((int(doping_time_sort[-1]),int(baking_time_sort[-1]+1)))
Z_top=np.zeros((int(doping_time_sort[-1]),int(baking_time_sort[-1]+1)))
Z_bottom=np.zeros((int(doping_time_sort[-1]),int(baking_time_sort[-1]+1)))
for i in range(0,len(doping_time)):
    Z_best[int(doping_time[i])-1][int(baking_time[i])]=polish_depth_best[i]
    Z_top[int(doping_time[i])-1][int(baking_time[i])]=polish_depth_top[i]
    Z_bottom[int(doping_time[i])-1][int(baking_time[i])]=polish_depth_bottom[i]
    #print(int(doping_time[i]),int(baking_time[i]-1))
Z_best=np.array(Z_best)
Z_top=np.array(Z_top)
Z_bottom=np.array(Z_bottom)

fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
#以下注释掉的是scatter的画图代码
#ax.scatter(doping_time,baking_time,polish_depth_best,label="best")
#ax.scatter(doping_time,baking_time,polish_depth_top,label="top")
#ax.scatter(doping_time,baking_time,polish_depth_bottom,label="bottom")
#以下是surface的画图代码
#XX是
YY,XX=np.meshgrid(np.array(range(0,int(baking_time_sort[-1])+1)),np.array(range(1,int(doping_time_sort[-1])+1)))
ax.plot_surface(XX,Z_best,YY,alpha=0.9)
ax.plot_surface(XX,Z_top,YY,alpha=0.5)
ax.plot_surface(XX,Z_bottom,YY,alpha=0.5)

#ax.contour(XX,YY,Z_best,zdir='y', offset=-6,cmap="rainbow")

ax.set_xlabel("doping time, min")
ax.set_zlabel("baking time, min")
ax.set_ylabel("polishing depth, μm")
plt.title("800℃ 3.3Pa ")
ax.legend(loc=0,ncol=1)
plt.show()
'''
plt.figure(1)
plt.tick_params(labelsize=10)
plt.plot(doping_time,polish_depth_best,label="best")
plt.plot(doping_time,polish_depth_top,label="top")
plt.plot(doping_time,polish_depth_bottom,label="bottom")
plt.fill_between(doping_time, polish_depth_top,polish_depth_bottom,facecolor='green', alpha=0.3)
plt.xlabel("doping time, min")
plt.ylabel("polishing depth, μm")
plt.title("800℃ 3.3Pa "+baking_time+" baking")
plt.legend(loc=0,ncol=1)
plt.show()
'''
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