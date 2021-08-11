import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt


thickness_a=500*10**-6
D=6.8*10**-14
s=0.075
p=4.7
kb=1.38*10**-23
T=800+273
mass=28*10**-3/(6.02*10**23)
J=s*p*(1/(2*m.pi*kb*T*mass))**0.5*2


total_t=3600*3 #120 seconds
n_t=total_t*10  #120 seconds
n_x=1000
delta_x=thickness_a/n_x
delta_t=total_t/n_t
Cs=2*J*(D*120)**0.5/D
alpha=D*delta_t/(2*delta_x**2)
#u_max=0.093*6.02*10**23*10**6
u_max=1*10**28

print("%.6e"%J)
print("%.6e"%(J*D/delta_x))
print("%.6e"%alpha)




#Initialize u and x

x=[]
y=[]
u=[[]]
for j in range(n_x):
  u[0].append(0)
  x.append((j+1)*delta_x)
  y.append(Cs*m.erfc((j+1)*delta_x/(2*(D*120)**0.5)))

#Set matrix A
A=[]
for i in range(n_x):
  A_row=[]
  for j in range(n_x):
    if(i==0):
      if(j==0):
        A_row.append(1+alpha)
      if(j==1):
        A_row.append(-alpha)
      if(j>1):
        A_row.append(0)
    else:
      if(j<i-1):
        A_row.append(0)
      if(j==i-1):
        A_row.append(-alpha)
      if(j==i):
        A_row.append(1+2*alpha)
      if(j==i+1):
        A_row.append(-alpha)
      if(j>i+1):
        A_row.append(0)
  #print(A_row)
  A.append(A_row)
A=np.matrix(A)

#Caculate the distridution of Nitrogen content
for k in range(n_t):
  print(k)
  alpha=D*delta_t/(2*delta_x**2)

  
  for i in range(n_x):
    if(u[k][i]<=2*u_max/5):
      alpha=1.37*10**(-13)*delta_t/(2*delta_x**2)
    
    for j in range(n_x):
      if(i==0):
        if(j==0):
          A[i,j]=1+alpha
        if(j==1):
          A[i,j]=-alpha
      else:
        if(j==i-1):
          A[i,j]=-alpha
        if(j==i):
          A[i,j]=1+2*alpha
        if(j==i+1):
          A[i,j]=-alpha
    #print(A_row)
  A=np.matrix(A)

  
  # Set matrix B
  B=[]
  for i in range(n_x):
    alpha=D*delta_t/(2*delta_x**2)
    if(u[k][i]<=2*u_max/5):
      alpha=1.37*10**(-13)*delta_t/(2*delta_x**2)
    J=s*p*(1/(2*m.pi*kb*T*mass))**0.5*2
    if(i==0):
      J=J*delta_x*(u_max-u[k][0])/(delta_x*u_max+D*J)
      B.append([u[k][0]*(1-alpha)+u[k][1]*alpha+2*alpha*J*delta_x/D])
    else:
      if(i==n_x-1):
        B.append([alpha*u[k][i-1]+(1-2*alpha)*u[k][i]])
      else:
        B.append([alpha*u[k][i-1]+(1-2*alpha)*u[k][i]+alpha*u[k][i+1]])
  B=np.matrix(B)
  result=np.linalg.solve(A,B)
  u_n=[]
  for i in range(n_x):
    u_n.append(result[i,0])
  u.append(u_n)
for i in range(10):

  print("%.6e"%u[n_t][i],end=" ")
print(" ")

#Caculate 6min without Nitrogen

u_6min=[]
u_6min.append(u[n_t])

total_t=1 #120 seconds
n_t=total_t*1  
n_x=1000
delta_x=thickness_a/n_x
delta_t=total_t/n_t
alpha=D*delta_t/(2*delta_x**2)
J=0

for k in range(n_t):
  
  
  alpha=D*delta_t/(2*delta_x**2)

  A=[]
  for i in range(n_x):
    if(u_6min[k][i]<=2*u_max/5):
      alpha=1.37*10**(-13)*delta_t/(2*delta_x**2)
    A_row=[]
    for j in range(n_x):
      if(i==0):
        if(j==0):
          A_row.append(1+alpha)
        if(j==1):
          A_row.append(-alpha)
        if(j>1):
          A_row.append(0)
      else:
        if(j<i-1):
          A_row.append(0)
        if(j==i-1):
          A_row.append(-alpha)
        if(j==i):
          A_row.append(1+2*alpha)
        if(j==i+1):
          A_row.append(-alpha)
        if(j>i+1):
          A_row.append(0)
    #print(A_row)
    A.append(A_row)
  A=np.matrix(A)


  # Set matrix B
  B=[]
  for i in range(n_x):
    if(i==0):
      B.append([u_6min[k][0]*(1-alpha)+u_6min[k][1]*alpha+2*alpha*J*delta_x/D])
    else:
      if(i==n_x-1):
        B.append([alpha*u_6min[k][i-1]+(1-2*alpha)*u_6min[k][i]])
      else:
        B.append([alpha*u_6min[k][i-1]+(1-2*alpha)*u_6min[k][i]+alpha*u_6min[k][i+1]])
  B=np.matrix(B)
  result=np.linalg.solve(A,B)
  u_n=[]
  for i in range(n_x):
    u_n.append(result[i,0])
  u_6min.append(u_n)
for i in range(10):

  print("%.6e"%u_6min[n_t][i],end=" ")
print(" ")


#plt.axis([0,30*10**-6,10**21,10**29])
plt.xlim(0,300*10**-6)
plt.semilogy(x,u[3600*3*10])
plt.semilogy(x,u_6min[n_t])
#plt.plot(x,y)
plt.show()
  



      