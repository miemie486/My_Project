import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy.fftpack import ifft
import csv

#读取文件
#filename="9cell1"
#f1=open(filename+".txt")

#f1=open("20231116 5 cell 650MHz_phase_pi mode_32800_6.txt")
f1=open("20231117_Fourier_filter_data.csv")

data_f1=f1.readlines()
length=len(data_f1)
dataset=[]
for i in range(length):
    data_tuple=(float(data_f1[i].split(',')[0]),float(data_f1[i].split(',')[1]))
    dataset.append(data_tuple)

dataset=list(set(dataset))
dataset=sorted(dataset,key=lambda x: x[0])
length=len(dataset)

f_out=open('20231117.csv','w', newline='')
writer=csv.writer(f_out)

data_x=[]
data_y=[]

for i in range(length):
    data_x.append(dataset[i][0])
    data_y.append(dataset[i][1])
    writer.writerow([dataset[i][0],dataset[i][1]])

f_out.close()

#求极值

data_valley_y=[]
data_valley_x=[]
cache_parameter=100
for i in range(cache_parameter,length-cache_parameter):
    cache=data_y[i-cache_parameter:i+cache_parameter]
    if((data_y[i-1]-data_y[i])*(data_y[i]-data_y[i+1])<0):
        cache_=sorted(cache)
        if(cache_[0]>=data_y[i]):
            data_valley_y.append(data_y[i])
            data_valley_x.append(data_x[i])
            #print(cache,data_y[i])

'''
plt.figure(1)
plt.scatter(data_x,data_y,color="red")
plt.plot(data_x,data_y,label="Flatten")
'''

print("background right=",np.sum(data_y[-1000:-1])/1000)
print("background left=",np.sum(data_y[0:1000])/1000)
for i in range(len(data_valley_x)):
    print('极值点:(',data_valley_x[i],',',data_valley_y[i],')')

plt.figure(3)
plt.scatter(data_valley_x,data_valley_y,color="red")
plt.plot(data_x,data_y)
#plt.plot(data_x,data_y)


plt.show()
