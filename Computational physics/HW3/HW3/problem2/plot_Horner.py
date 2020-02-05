import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

a=0.3
f1=open("dataHorner.txt")
f2=open("datafloat1.txt")
f3=open("datafloat2.txt")
plt.figure(1)
t111=f1.readlines()          #从这里是读取txt文件数据的开始
p=len(t111)                   #这一段代码处理的是一个有两列数据的txt文件
i=0                           #将txt文件里面的每一列储存到一个数组当中
t211x=[]
t211y=[]
while i<=p-1:
    t311=t111[i].split()
    t511x=float(t311[0])
    t211x.append(t511x)
    t511y=float(t311[1])      #这里是结束
    
    t211y.append(t511y)
    i=i+1
t211x=np.array(t211x)
t211y=np.array(t211y)
plt.scatter(t211x,t211y)
plt.xlim(1-a,1+a)
plt.ylim(-10**-5,10**-5)
plt.title("dataHorner")

plt.figure(2)
t112=f2.readlines()          #读取文件
p=len(t112)
i=0
t212x=[]
t212y=[]
while i<=p-1:
    t312=t112[i].split()
    t512x=float(t312[0])
    t212x.append(t512x)
    t512y=float(t312[1])
    
    t212y.append(t512y)
    i=i+1
t212x=np.array(t212x)
t212y=np.array(t212y)
plt.scatter(t212x,t212y)
plt.xlim(1-a,1+a)
plt.ylim(-10**-6,10**-6)
plt.title("datafloat1")

plt.figure(3)
t121=f3.readlines()          #读取文件
p=len(t121)
i=0
t221x=[]
t221y=[]
while i<=p-1:
    t321=t121[i].split()
    t521x=float(t321[0])
    t221x.append(t521x)
    t521y=float(t321[1])
    
    t221y.append(t521y)
    i=i+1
t221x=np.array(t221x)
t221y=np.array(t221y)
plt.scatter(t221x,t221y)
plt.xlim(1-a,1+a)
plt.ylim(-10**-6,10**-6)
plt.title("datafloat2")

plt.show()

f1.close()
f2.close()
f3.close()