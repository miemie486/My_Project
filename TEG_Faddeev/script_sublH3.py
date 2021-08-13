import os
import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
from numba import jit

if not os.path.exists('jobs/data'):
    os.makedirs('jobs/data')
data_final=[]
for j in [-0.02,-0.01,0.01,0.02]:

    command1='cd jobs&&../bin/sublH3 '
    command2=str(0.01*j)
    command3=' > data/test.txt'
    command=command1+command2+command3
    #print(command)
    os.system(command)
    f=open('jobs/data/test.txt')
    data0=f.readlines()
    data1=[]
    data=[]  #data才是最终提取出来的数据

    l=0   #这个是中间变量，用于提取数据
    for string in data0:
        if(l==1):
            data1.append(string)
        if(string.find('Mesh Points         Binding Energy (MeV)')!=-1):
            l=1

    for i in data1:  #t0,t1是仅仅是这个循环的中间变量
        t0=i.split()
        t1=list(map(float,t0))
        data.append(t1)
    data_final.append(data[0][1])

#print(data_final)

        


