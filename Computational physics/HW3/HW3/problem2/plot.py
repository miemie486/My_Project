import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt


f12=open("datafloat2.txt")
f21=open("datadouble1.txt")
f22=open("datadouble2.txt")
f31=open("dataquad1.txt")
f32=open("dataquad2.txt")
a=0.3      #点集的区间为1+-a
b=0.04
c=0.0004
plt.figure(1)
ax11=plt.subplot(2,3,1)
ax12=plt.subplot(2,3,4)
ax21=plt.subplot(2,3,2)
ax22=plt.subplot(2,3,5)
ax31=plt.subplot(2,3,3)
ax32=plt.subplot(2,3,6)
#floatdata1画图开始
f11=open("datafloat1.txt")
t111=f11.readlines()          #从这里是读取txt文件数据的开始
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
plt.sca(ax11)
plt.scatter(t211x,t211y)
plt.xlim(1-a,1+a)
plt.ylim(-10**-6,10**-6)
plt.title("datafloat1")
#floatdata1画图结束
#floatdata2画图开始
t112=f12.readlines()          #读取文件
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
plt.sca(ax12)
plt.scatter(t212x,t212y)
plt.xlim(1-a,1+a)
plt.ylim(-10**-6,10**-6)
plt.title("datafloat2")
#floatdata2画图结束
#doubledata1画图开始
t121=f21.readlines()          #读取文件
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
plt.sca(ax21)
plt.scatter(t221x,t221y)
plt.xlim(1-b,1+b)
plt.ylim(-10**-13,10**-13)
plt.title("datadouble1")
#doubledata1画图结束
#doubledata2画图开始
t122=f22.readlines()          #读取文件
p=len(t122)
i=0
t222x=[]
t222y=[]
while i<=p-1:
    t322=t122[i].split()
    t522x=float(t322[0])
    t222x.append(t522x)
    t522y=float(t322[1])
    
    t222y.append(t522y)
    i=i+1
t222x=np.array(t222x)
t222y=np.array(t222y)
plt.sca(ax22)
plt.scatter(t222x,t222y)
plt.xlim(1-b,1+b)
plt.ylim(-10**-13,10**-13)
plt.title("datadouble2")
#doubledata2画图结束
#quaddata1画图开始
t131=f31.readlines()          #读取文件
p=len(t131)
i=0
t231x=[]
t231y=[]
while i<=p-1:
    t331=t131[i].split()
    t531x=float(t331[0])
    t231x.append(t531x)
    t531y=float(t331[1])
    
    t231y.append(t531y)
    i=i+1
t231x=np.array(t231x)
t231y=np.array(t231y)
plt.sca(ax31)
plt.scatter(t231x,t231y)
plt.xlim(1-c,1+c)
plt.ylim(-10**-31,10**-31)
plt.title("dataquadmath1")
#quaddata1画图结束
#quaddata2画图开始
t132=f32.readlines()          #读取文件
p=len(t132)
i=0
t232x=[]
t232y=[]
while i<=p-1:
    t332=t132[i].split()
    t532x=float(t332[0])
    t232x.append(t532x)
    t532y=float(t332[1])
    
    t232y.append(t532y)
    i=i+1
t232x=np.array(t232x)
t232y=np.array(t232y)
plt.sca(ax32)
plt.scatter(t232x,t232y)
plt.xlim(1-c,1+c)
plt.ylim(-10**-31,10**-31)
plt.title("dataquadmath2")
#quaddata2画图结束
plt.show()
f11.close()
f12.close()
f21.close()
f22.close()
f31.close()
f32.close()

