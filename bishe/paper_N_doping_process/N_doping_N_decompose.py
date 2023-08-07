import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from numba import jit
from numpy.core import umath



thickness_a=50*10**-6
s=0.075*0.00075*0.2
#s=0.0001*0.075*7.5
kb=1.38*10**-23
mass=28*10**-3/(6.02*10**23) #氮分子质量


concentration_ratio=3.27#4.41    #化学反应速率过渡区开始与结束点的浓度之比


#靠近表面那一层，在对数坐标轴的图上的斜率，主要是由化学反应速率常数影响的
#化学反应速率越小，斜率越小
#如果需要设定比较小的u_max或者较大的n_x，需要同时设定比较小的delta_t
#要不然在左边界会出现奇异点
'''
parameters=[T,p,C_N,u_max_1,umax_2,N_x,N_t,total_t,D]

T:温度      默认800+273单位K
p:氮气压强      默认3.3单位Pa
C_N:氮气化学反应速率        没有默认，随便乱改
u_max_1:开始生成氮化物的氮浓度      默认2.79e25
u_max_2:氮的最高浓度        默认2.46e28 因为Nb2N的密度为8.02g/cm^3,就能算出此时的N浓度
N_x:长度上的格子数，也是矩阵大小        默认200 
N_t:时间上的步数        默认total_t*10
total_t:总时间      看掺杂工艺，默认2*60_6*60工艺
D:氮在铌里的扩散系数，要随着T的改变一起变   默认800+273K情况下D=6.8*10**-14


'''

@jit
def set_matrix_A(parameters,u,u_total):        #构建矩阵A
    T=parameters[0]
    p=parameters[1]
    C_N=parameters[2]
    u_max_1=parameters[3]
    u_max_2=parameters[4]
    N_x=int(parameters[5])
    N_t=int(parameters[6])
    total_t=parameters[7]
    D=parameters[8]

    delta_x=thickness_a/N_x
    delta_t=total_t/N_t
    alpha=D*delta_t/(2*delta_x**2)
    
    #构建矩阵A
    i=0
    A=np.zeros(shape=(N_x,N_x))
    total_u=np.zeros(N_x)
    while(i<N_x):

        #对化学反应系数随浓度变化的修正
        if(u[i]<=u_max_1 or u[i]>=u_max_2 or total_u[i]>=u_max_2):
            C_N=0
        elif(u[i]>=u_max_1 and u[i]<=concentration_ratio*u_max_1 ):
            C_N=parameters[2]*m.log(u[i]/u_max_1)/m.log(concentration_ratio)
            #C_N=0
        elif(u[i]<=u_max_2 and u[i]>=0.1*u_max_2):
            C_N=parameters[2]*m.log10(u_max_2/u[i])
            #C_N=0
        else:
            C_N=parameters[2]

        if(i==0):
            A[i,0]=1+C_N*delta_t/2+alpha
            A[i,1]=-alpha
        elif(i<N_x-1):
            A[i,i-1]=-alpha
            A[i,i]=1+2*alpha+C_N*delta_t/2
            A[i,i+1]=-alpha
        else:
            A[i,i-1]=-alpha
            A[i,i]=1+2*alpha+C_N*delta_t/2
        
        total_u[i]=C_N*(u[i]-u_max_1)*delta_t
        i=i+1
    return A  

@jit
def set_matrix_B_and_reaction_u(parameters,u,_total_u):
    T=parameters[0]
    p=parameters[1]
    C_N=parameters[2]
    u_max_1=parameters[3]
    u_max_2=parameters[4]
    N_x=int(parameters[5])
    N_t=int(parameters[6])
    total_t=parameters[7]
    D=parameters[8]

    J_0=2*s*p*(1/(2*m.pi*kb*T*mass))**0.5  #每秒吸附在表面的氮浓度
    delta_x=thickness_a/N_x
    delta_t=total_t/N_t
    alpha=D*delta_t/(2*delta_x**2)

    #构建矩阵B
    i=0
    B=np.zeros(N_x)
    reaction_u=np.zeros(N_x)
    while(i<N_x):
        '''
        if(u[i]<=u_max_1 or u[i]>=u_max_2):
            C_N=0
        elif(u[i]>=u_max_1 and u[i]<=5*u_max_1):
            C_N=parameters[2]*(u[i]-u_max_1)/u_max_1
            #C_N=0
        elif(u[i]<=u_max_2 and u[i]>=0.2*u_max_2):
            C_N=parameters[2]*(u_max_2-u[i])/(0.8*u_max_2)
            #C_N=0
        else:
            C_N=parameters[2]
        '''
        #对化学反应系数随浓度变化的修正
        if(u[i]<=u_max_1 or u[i]>=u_max_2 or _total_u[i]>=u_max_2):
            C_N=0
        elif(u[i]>=u_max_1 and u[i]<=concentration_ratio*u_max_1):
            C_N=parameters[2]*m.log(u[i]/u_max_1)/m.log(concentration_ratio)
            #C_N=0
        elif(u[i]<=u_max_2 and u[i]>=0.1*u_max_2):
            C_N=parameters[2]*m.log10(u_max_2/u[i])
            #C_N=0
        else:
            C_N=parameters[2]


        #test如果totalu>umax1，u[i]固定
        if(_total_u[i]>u_max_1 and p==0):
            u[i]=u_max_1*0.9
        #test如果totalu>umax1，u[i]固定


        if(i==0):
            J=J_0*(1-_total_u[0]/u_max_2)#(u_max_2*D-u[0]*D)/(u_max_2*D+J_0*delta_x)
            B[i]=2*alpha*J*delta_x/D+(1-C_N*delta_t/2-alpha)*u[i]+alpha*u[i+1]+C_N*u_max_1*delta_t
        elif(i<N_x-1):
            B[i]=alpha*u[i-1]+(1-C_N*delta_t/2-2*alpha)*u[i]+alpha*u[i+1]+C_N*u_max_1*delta_t
        else:
            B[i]=alpha*u[i-1]+(1-C_N*delta_t/2-2*alpha)*u[i]+C_N*u_max_1*delta_t
        
        reaction_u[i]=C_N*(u[i]-u_max_1)*delta_t
        if(u[i]<=u_max_1 and _total_u[i]>0):
            B[i]=B[i]
            #B[i]=B[i]+C_N*(u_max_1-u[i])*delta_t
            #reaction_u[i]=-1.03*(u_max_1-u[i])


        


        i=i+1

    return [B,reaction_u]


def solve_equ(parameters,u,_total_u):

    T=parameters[0]
    C_N=parameters[2]
    N_x=int(parameters[5])
    N_t=int(parameters[6])
    total_t=parameters[7]
    D=parameters[8]


    i=0
    total_u=_total_u
    while(i<N_t):
        A=set_matrix_A(parameters,u,total_u)
        B_and_reaction_u=set_matrix_B_and_reaction_u(parameters,u,total_u)
        B=B_and_reaction_u[0]
        result_u=np.linalg.solve(A,B)
        u=result_u
        total_u=total_u+B_and_reaction_u[1]
        print(i*total_t/N_t,u[0],D*(u[0]-u[1])/(thickness_a/N_x),total_u[0])
        #print(B_and_reaction_u[0][0])
        i=i+1

    return [u,total_u]

if __name__ == "__main__":
    N_x=400
    total_t=20*60
    N_t=total_t*10
    C_N=0.4
    parameters=[800+273,3.3,C_N,2.79e25,2.85e27,N_x,N_t,total_t,6.8*10**-14]
    parameters=np.array(parameters)
    u=np.zeros(N_x)

    start = timer()

    parameters=[800+273,3.3,1.5,2.79e25,2.85e27,500,100,1,6.8*10**-14]
    parameters=np.array(parameters)

    result=solve_equ(parameters,u)

    parameters=[800+273,3.3,1.5,2.79e25,2.85e27,500,1200*5,1200,6.8*10**-14]
    parameters=np.array(parameters)

    result=solve_equ(parameters,result[0])

    X=np.linspace(0,thickness_a,500)
    plt.figure(1)
    plt.xlim(0,5*10**-6)
    plt.ylim(1e23,1e28)
    plt.semilogy(X,result[0]+1e18)
    #plt.semilogy(X,result[1])
    print("_____________________")
    print("计算耗时：",timer()-start)
    plt.show()
