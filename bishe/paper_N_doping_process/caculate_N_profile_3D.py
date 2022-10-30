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

if __name__ == "__main__":
    start = timer()
    data=importfile('data_800_20min.txt')

    C_N=3.48
    N_x=400

    X=np.linspace(0,N_D.thickness_a,N_x)
    np.save("data\X",X)

    #单位min
    #doping_and_baking_time=[[2,0],[2,6],[3,60],[10,0],[20,0],[20,30]]
    #doping_and_baking_time=[[0.2,0],[0.4,0],[0.6,0],[0.8,0],[1.5,0],[2.5,0],[3.5,0],[4.5,0]]
    '''
    doping_and_baking_time=[]
    for i in range(1,31):
        doping_and_baking_time.append([i,0])
    print(doping_and_baking_time)
    '''
    doping_time_range=range(1,31)
    baking_time_range=range(1,31)

    doping_but_no_baking_result=np.zeros(N_x)

    for i in doping_time_range:
        start_local = timer()
        #doping 过程
        doping_time=60
        N_doping_time=doping_time*10

        parameters=[800+273,3.3,C_N,2.79e25,2.46e28,N_x,N_doping_time,doping_time,6.8*10**-14]
        parameters=np.array(parameters)

        u=doping_but_no_baking_result
        result_d=N_D.solve_equ(parameters,u)

        doping_but_no_baking_result=result_d[0]

        name="data_3D\data_800_"+str(i)+"N0min_theory"
        np.save(name,result_d[0])

        doping_and_baking_result=doping_but_no_baking_result
        #baking 过程
        for j in baking_time_range:
            start_local2 = timer()
            baking_time=60
            N_baking_time=baking_time*10

            parameters=[800+273,0,C_N,2.79e25,2.46e28,N_x,N_baking_time,baking_time,6.8*10**-14]
            parameters=np.array(parameters)

            result_d_b=N_D.solve_equ(parameters,doping_and_baking_result)
            doping_and_baking_result=result_d_b[0]

            name="data_3D\data_800_"+str(i)+"N"+str(j)+"min_theory"
            
            np.save(name,result_d_b[0])
            print("_____________________")
            print(name+"计算耗时：",timer()-start_local2)
        print("_____________________")
        print(name+"计算耗时：",timer()-start_local)

    print("_____________________")
    print("总计算耗时：",timer()-start)