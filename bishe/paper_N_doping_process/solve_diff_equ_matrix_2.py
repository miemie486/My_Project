import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt


thickness_a=50*10**-6
D=6.8*10**-14
s=0.075
p=3.3
kb=1.38*10**-23
T=800+273
mass=28*10**-3/(6.02*10**23)
J=s*p*(1/(2*m.pi*kb*T*mass))**0.5*2

#靠近表面那一层，在对数坐标轴的图上的斜率，主要是由化学反应速率常数影响的
#化学反应速率越小，斜率越小

k_N_1=5   #氮的反应速率常数
k_N=5

#如果需要设定比较小的u_max或者较大的n_x，需要同时设定比较小的delta_t
#要不然在左边界会出现奇异点

total_t=2*60 #120 seconds
n_t=int(total_t*10)  #120 seconds
n_x=200
delta_x=thickness_a/n_x
delta_t=total_t/n_t
Cs=2*J*(D*120)**0.5/D
alpha=D*delta_t/(2*delta_x**2)
#u_max=0.093*6.02*10**23*10**6
u_max=0.15*1e28
#u_max=1

print("%.6e"%J)
print("%.6e"%(J*D/delta_x))
print("%.6e"%alpha)

name1="N_2min_3.3.txt"
name2="N_2min_6min_3.3.txt"




#Initialize u and x

x=[]
y=[]
u=[[]]
u_total_N=[]
for j in range(n_x):
  u[0].append(0)
  u_total_N.append(0)
  x.append((j+1)*delta_x)
  y.append(Cs*m.erfc((j+1)*delta_x/(2*(D*120)**0.5)))

#Set matrix A
A=[]
for i in range(n_x):
  A_row=[]
  for j in range(n_x):
    if(i==0):
      if(j==0):
        A_row.append(1+alpha+k_N*delta_t/2)
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
        A_row.append(1+2*alpha+k_N*delta_t/2)
      if(j==i+1):
        A_row.append(-alpha)
      if(j>i+1):
        A_row.append(0)
  #print(A_row)
  A.append(A_row)
A=np.matrix(A)

#Caculate the distridution of Nitrogen content
for k in range(n_t):
  #rewrite A
  alpha=D*delta_t/(2*delta_x**2)
  for i in range(n_x):
    """
    if(u[k][i]<=2*u_max/50):
      #alpha=1.37*10**(-13)*delta_t/(2*delta_x**2)
      alpha=D*delta_t/(2*delta_x**2)
    """

    if(u[k][i]<=0.01*u_max or u_total_N[i]>=5e28):
      k_N=0
      D=6.8*10**-14
    else:
      
      #print(i)

      k_N=k_N_1
      D=6.8*10**-14

    if(i==0):
      A[i,0]=1+alpha+k_N*delta_t/2
      A[i,1]=-alpha
    elif(i<n_x-1):
      A[i,i-1]=-alpha
      A[i,i]=1+2*alpha+k_N*delta_t/2
      A[i,i+1]=-alpha
    else:
      A[i,i-1]=-alpha
      A[i,i]=1+2*alpha+k_N*delta_t/2
    #print(A_row)
  A=np.array(A)

  
  # Set matrix B
  B=np.empty((n_x,1))
  for i in range(n_x):
    alpha=D*delta_t/(2*delta_x**2)
    """
    if(u[k][i]<=2*u_max/50):
      #alpha=1.37*10**(-13)*delta_t/(2*delta_x**2)
      alpha=D*delta_t/(2*delta_x**2)
    """


    if(u[k][i]<=0.01*u_max or u_total_N[i]>=5e28):
      k_N=0
      D=6.8*10**-14
    else:
      k_N=k_N_1
      D=6.8*10**-14/5



    J=s*p*(1/(2*m.pi*kb*T*mass))**0.5*2
    if(i==0):
      J=J*delta_x*(u_max-u[k][0])/(delta_x*u_max+D*J)

      print(k/10,"%.6e"%J,"%.6e"%u[k][0])

      B[i,0]=u[k][0]*(1-alpha)+u[k][1]*alpha+2*alpha*J*delta_x/D
    else:
      if(i==n_x-1):
        B[i,0]=alpha*u[k][i-1]+(1-2*alpha-k_N*delta_t/2)*u[k][i]
      else:
        B[i,0]=alpha*u[k][i-1]+(1-2*alpha-k_N*delta_t/2)*u[k][i]+alpha*u[k][i+1]
  B=np.array(B)
  result=np.linalg.solve(A,B)
  u_n=[]

  for i in range(n_x):
    u_n.append(result[i,0])
    if(u[k][i]<=0.01*u_max or u_total_N[i]>=5e28):
      k_N=0
      u_total_N[i]=u_total_N[i]+k_N*result[i,0]*delta_t
    else:
      k_N=k_N_1
      u_total_N[i]=u_total_N[i]+k_N*result[i,0]*delta_t

  u.append(u_n)
for i in range(10):

  print("%.6e"%u[n_t][i],end=" ")
print(" ")

#设定最小值
for i in range(n_x):
  u_n.append(result[i,0])
  u_total_N[i]=u_total_N[i]+1e-4*u_max+result[i,0]




f1=open(name1,mode="w+")
for i in range(n_x):
  f1.write(str(x[i]))
  f1.write(" ")
  f1.write(str(u[n_t][i]))
  f1.write("\n")
f1.close()

#Caculate 6min without Nitrogen

u_6min=[]
u_6min.append(u[n_t])

total_t=1 #360 seconds
n_t=total_t*10



n_t=1



delta_x=thickness_a/n_x
delta_t=total_t/n_t
alpha=D*delta_t/(2*delta_x**2)
J=0

for k in range(n_t):
  
  
  #print(k/10)
  alpha=D*delta_t/(2*delta_x**2)
  for i in range(n_x):
    """
    if(u_6min[k][i]<=2*u_max/50):
      #alpha=1.37*10**(-13)*delta_t/(2*delta_x**2)
      alpha=D*delta_t/(2*delta_x**2)
    """

    if(i==0):
      A[i,0]=1+alpha
      A[i,1]=-alpha
    elif(i<n_x-1):
      A[i,i-1]=-alpha
      A[i,i]=1+2*alpha
      A[i,i+1]=-alpha
    else:
      A[i,i-1]=-alpha
      A[i,i]=1+2*alpha
    #print(A_row)
  A=np.matrix(A)


  # Set matrix B
  B=np.empty((n_x,1))
  for i in range(n_x):
    alpha=D*delta_t/(2*delta_x**2)
    """
    if(u_6min[k][i]<=2*u_max/50):
      #alpha=1.37*10**(-13)*delta_t/(2*delta_x**2)
      alpha=D*delta_t/(2*delta_x**2)
    """
    J=s*p*(1/(2*m.pi*kb*T*mass))**0.5*2
    if(i==0):
      #J=-J*delta_x*u_6min[k][0]/(delta_x*u_max+D*J)
      J=0
      print(k/10,"%.6e"%J)
      B[i,0]=u_6min[k][0]*(1-alpha)+u_6min[k][1]*alpha+2*alpha*J*delta_x/D
    else:
      if(i==n_x-1):
        B[i,0]=alpha*u_6min[k][i-1]+(1-2*alpha)*u_6min[k][i]
      else:
        B[i,0]=alpha*u_6min[k][i-1]+(1-2*alpha)*u_6min[k][i]+alpha*u_6min[k][i+1]
  B=np.matrix(B)
  result=np.linalg.solve(A,B)
  u_n=[]
  for i in range(n_x):
    u_n.append(result[i,0])

  u_6min.append(u_n)
for i in range(10):

  print("%.6e"%u_6min[n_t][i],end=" ")
print(" ")
f2=open(name2,mode="w+")
for i in range(n_x):
  f2.write(str(x[i]))
  f2.write(" ")
  f2.write(str(u_6min[n_t][i]))
  f2.write("\n")
f2.close()

#plt.axis([0,30*10**-6,10**21,10**29])
plt.figure(1)
plt.xlim(0,30*10**-6)
#plt.plot(x,u[120*10])
plt.semilogy(x,u_6min[n_t])
plt.figure(2)
plt.xlim(0,10*10**-6)
plt.ylim(1e23,1e29)
plt.semilogy(x,u_6min[n_t],label="N alpha phase")
plt.semilogy(x,u_total_N,label="N beta phase")
plt.legend(loc=0,ncol=1)
#plt.plot(x,y)
plt.show()
  



      