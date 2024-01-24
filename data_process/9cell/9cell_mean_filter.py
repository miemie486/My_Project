import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.fftpack import fft


filename="9cell2"
f1=open(filename+".txt")
data_f1=f1.readlines()
length=len(data_f1)

data_x=[]
data_y=[]

for i in range(length):
    data_x.append(float(data_f1[i].split()[0]))
    data_y.append(float(data_f1[i].split()[1]))

n=80


n_left_array=[]
n_right_array=[]

for i in range(n):
    n_left_array.append(data_y[0])
    n_right_array.append(data_y[-1])

_data_y=np.array(n_left_array+data_y+n_right_array)
sum_array=[]
for i in range(length):
    sum_array.append(_data_y[i+n-n:i+n+n])


f_out=open(filename+"_mean_filter_data.txt",'w')

data_final_x=data_x
data_final_y=[]
error_2=0
for i in range(length):
    Flatten_y=np.sum(sum_array[i])/(2*n)
    data_final_y.append(Flatten_y)
    error_2=error_2+(_data_y[i]-Flatten_y)**2
    print(str(data_final_x[i]),str(Flatten_y),file=f_out)

f_out.close()
    
error=m.sqrt(error_2)
print("error=",error)

#以下注释掉的是查看频谱
'''
fft_y=abs(fft(data_y))
fft_y_final=abs(fft(data_final_y))
fft_x=np.linspace(1,p,p)
fft_error2=0
for i in range(len(fft_y)):
    fft_error2=fft_error2+(fft_y[i]-fft_y_final[i])**2
fft_error=m.sqrt(fft_error2)
print("fft error=",fft_error)
print("background right=",np.sum(data_y[-1000:-1])/1000)
print("background left=",np.sum(data_y[0:1000])/1000)

print("background right=",np.sum(data_final_y[-1000:-1])/1000)
print("background left=",np.sum(data_final_y[0:1000])/1000)
'''
#求极值
data_y_D=[]

for i in range(length):
    if(data_final_y[i]>0):
        data_y_D.append(0)
    else:
        data_y_D.append(data_final_y[i])

data_valley_y=[]
data_valley_x=[]
for i in range(length-2):
    if((data_final_y[i]-data_final_y[i+1])*(data_final_y[i+1]-data_final_y[i+2])<0):
        if(data_final_y[i+1]<0):
            data_valley_y.append(data_final_y[i+1])
            data_valley_x.append(data_final_x[i+1])
            print(data_final_x[i+1],data_final_y[i+1])

    #Dy_Dx.append((data_final_y[i]-data_final_y[i+1])/(data_final_x[i]-data_final_x[i+1]))
    #D_x.append(data_final_x[i])



plt.figure(1)
plt.scatter(data_x,data_y,color="red")
plt.plot(data_final_x,data_final_y,label="Flatten")

plt.figure(3)
plt.scatter(data_valley_x,data_valley_y,color="red")
plt.plot(data_final_x,data_y_D)
'''
plt.figure(2)
plt.ylim(0,1000)
plt.plot(fft_x,fft_y,color="red")
plt.plot(fft_x,fft_y_final)
'''


plt.show()
