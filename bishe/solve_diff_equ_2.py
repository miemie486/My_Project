import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt


thickness_a=100*10**-6
D=6.8*10**-14
s=0.075
p=4.7
kb=1.38*10**-23
T=800+273
mass=28*10**-3/(6.02*10**23)
J=s*p*(1/(2*m.pi*kb*T*mass))**0.5*2

total_t=120*3 #120 seconds
n_t=total_t*10  #120 seconds
n_x=30 
delta_x=thickness_a/n_x
delta_t=total_t/n_t

u=[[2*J*(D*120)**0.5/D*0.9]]
#print("%.6e"%0,"%.6e"%0,"%.6e"%(u[0][0]))
x=[0]
y=[2*J*(D*120)**0.5/D]
#Initialize u and x
for j in range(n_x):
  u[0].append(y[0]*m.erfc((j+1)*delta_x/(2*(D*120)**0.5)))
  x.append((j+1)*delta_x)
  y.append(y[0]*m.erfc((j+1)*delta_x/(2*(D*120)**0.5)))
  #print("%.6e"%0,"%.6e"%((j+1)*delta_x),"%.6e"%(u[0][j+1]))

def f():
  sum_=0
  for i in range(int(n_t)):
    u_time_fixed=[]
    
    sum=0
    for j in range(n_x):
      if(j==n_x-1):
        a1=(D*delta_t/delta_x**2)*(u[i][j]+0)+(1-(D*delta_t/delta_x**2))*u[i][j+1]
      else:
        a1=(D*delta_t/delta_x**2)*(u[i][j]+u[i][j+2])+(1-(D*delta_t/delta_x**2))*u[i][j+1]

      if(j==0):
        u_time_fixed.append(a1*0.9)
        sum=sum+(a1)*delta_x
        #print("%.6e"%((i+1)*delta_t),"%.6e"%((j)*delta_x),"%.6e"%(u_time_fixed[j]),"%.6e"%(u_time_fixed[j]-u[i][j]))
        u_time_fixed.append(a1)
      else:
        u_time_fixed.append(a1)
      sum=sum+(a1)*delta_x
      #print("%.6e"%((i+1)*delta_t),"%.6e"%((j+1)*delta_x),"%.6e"%(u_time_fixed[j+1]),"%.6e"%(u_time_fixed[j]-u[i][j]))
    u.append(u_time_fixed)
    print(i/10,"%.6e"%((sum-sum_)/(J*delta_t)),"%.6e"%sum)
    sum_=sum
  plt.plot(x,u[len(u)-1])
  plt.plot(x,y)
  plt.show()


f()
    






