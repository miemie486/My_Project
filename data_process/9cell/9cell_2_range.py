import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.fftpack import fft

f1=open("9cell.txt")
data_f1=f1.readlines()
p=len(data_f1)

data_x=[]
data_y=[]

for i in range(p):
    data_x.append(float(data_f1[i].split()[0]))
    data_y.append(float(data_f1[i].split()[1]))

n=100
for kk in range(300):
    n=kk

    n_left_array=[]
    n_right_array=[]

    for i in range(n):
        n_left_array.append(data_y[0])
        n_right_array.append(data_y[-1])

    _data_y=np.array(n_left_array+data_y+n_right_array)
    sum_array=[]
    for i in range(p):
        sum_array.append(_data_y[i+n-n:i+n+n])

    data_final_x=data_x
    data_final_y=[]
    error_2=0
    for i in range(p):
        Flatten_y=np.sum(sum_array[i])/(2*n)
        data_final_y.append(Flatten_y)
        error_2=error_2+(_data_y[i]-Flatten_y)**2
        
    error=m.sqrt(error_2)
 



    fft_y=abs(fft(data_y))
    fft_y_final=abs(fft(data_final_y))
    fft_x=np.linspace(1,p,p)
    fft_error2=0
    for i in range(len(fft_y)):
        fft_error2=fft_error2+(fft_y[i]-fft_y_final[i])**2
    fft_error=m.sqrt(fft_error2)

    print("error=",error,"fft error=",fft_error)



'''
plt.figure(1)
plt.scatter(data_x,data_y,color="red")
plt.plot(data_final_x,data_final_y,label="Flatten")

plt.figure(2)
plt.ylim(0,1000)
plt.plot(fft_x,fft_y,color="red")
plt.plot(fft_x,fft_y_final)

'''

plt.show()
