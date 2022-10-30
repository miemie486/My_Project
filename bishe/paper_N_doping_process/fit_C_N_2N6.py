import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from numba import jit
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


#对于计算氮分布的具体过程得在这里调控

def sum_distence_square(data_exp_x,data_exp_y,C_N):

    N_x=400
    total_t=20*60
    N_t=total_t*10

    parameters=[800+273,3.3,C_N,2.79e25,2.765e28,N_x,N_t,total_t,6.8*10**-14]
    parameters=np.array(parameters)
    u=np.zeros(N_x)
    result=N_D.solve_equ(parameters,u)

    data_x=np.linspace(0,N_D.thickness_a,N_x)
    data_y=result[0]

    a=0
    i=0
    while(i<len(data_exp_x)):
        a=a+(data_exp_y[i]-np.interp(data_exp_x[i],data_x,data_y))**2
        i=i+1
    return a


    
if __name__ == "__main__":
    start = timer()
    data=importfile('data_800_20min.txt')

    C_N=3.48
    N_x=400
    total_t=2*60
    N_t=total_t*100

    parameters=[800+273,3.5,C_N,2.79e25,2.46e28,N_x,N_t,total_t,6.8*10**-14]
    parameters=np.array(parameters)

    u=np.zeros(N_x)
    result=N_D.solve_equ(parameters,u)

    '''
    total_t=60*60
    N_t=total_t*10

    parameters=[800+273,1e-10,C_N,2.79e25,2.46e28,N_x,N_t,total_t,6.8*10**-14]
    parameters=np.array(parameters)

    result=N_D.solve_equ(parameters,result[0])
    '''
    X=np.linspace(0,50*10**-6,400)
    
    '''
    a=0
    i=0
    while(i<42):
        a=a+(m.log10(data[1][i])-m.log10(np.interp(data[0][i],X,result[0])))**2
        i=i+1
    print(a)
    '''
    np.save("_test_data_800_10N0min_3.5Pa_theory",result[0])
    plt.figure(1)
    plt.xlim(0,5*10**-6)
    plt.ylim(1e23,2.46e28)
    plt.semilogy(X,result[0]+1e18,label="Theory")
    plt.semilogy(data[0],data[1],label="Experiment")
    plt.legend(loc=0,ncol=1)
    print("计算耗时：",timer()-start)
    plt.show()
    

    #print(sum_distence_square(data[0],data[1],1.5))