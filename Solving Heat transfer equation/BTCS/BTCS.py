#===================================
#                BTCS
#===================================
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

m=40                          # number of horizontal divisions     
n=40                          # number of vertical divisions  

l=1.                          # length of the domain
h=1.                          # height of the domain
dx=l/m
dy=h/n
be=dy/dx
dt=5                          # time step
Tf=30000                      # final time
alpha=23.1*1e-6               # diffusivity of the iron
S_T=5e-3                      # heat source
w=1.1                         # over relaxation factor
cr1=1e-3                      # convergence criteria for reaching steady state solution
cr2=1e-6                      # convergence criteria for Gauss_Seidel iteration


# initial condition
T=np.zeros((m+1,n+1));T_old=np.zeros((m+1,n+1));TT=np.zeros((m+1,n+1));

# boundary conditions
for i in range(m+1):
    T[i,0]=20
    T[i,n]=40

for j in range(n+1):
    T[m,j]=10



#=============================== solution ==============================
t=0
err1=1
while err1>cr1 and t<Tf:
    t=t+dt
    
    # updating boundary conditions
    for j in range(n+1):
        T[0,j]=T[1,j]
        
    T_old[:,:]=T[:,:]
    
    
    # BTCS loop  
    k=0
    err2=1
    while err2>cr2:
        k=k+1
        
        # updating boundary conditions
        for j in range(n+1):
            T[0,j]=T[1,j]
        
        TT[:,:]=T[:,:]
        
        err2=0
        for i in range(1,m):
            for j in range(1,n):
                T[i,j]=(1-w)*T[i,j]+(w/ (2*(1+be**2)+dy**2/dt/alpha) )*( T[i,j-1]+T[i,j+1] +be**2*( T[i-1,j]+T[i+1,j] )  +S_T/alpha*dy**2  +T_old[i,j]*dy**2/dt/alpha )
                err2=err2+math.sqrt( (T[i,j]-TT[i,j])**2/((m-1)*(n-1)) )        
     
        
    err1=0
    for i in range(1,m):
        for j in range(1,n):
            err1=err1+math.sqrt( (T[i,j]-T_old[i,j])**2/((m-1)*(n-1)) )    
            
#    plt.scatter(t,err1)    
#    plt.xlabel('time(s)');plt.ylabel('norm of the error')
#    plt.title('convergence history')
#    plt.yscale("log")
#    plt.pause(0.0001)
    
    # printing time and second norm of the error
    print('{:10.2f}{:15.2e}{:6d}'.format(t,err1,k))
    
# updating boundary conditions
for j in range(n+1):
    T[0,j]=T[1,j]
        
        

#=============================== results ===============================
# printing temperature at the middle of the plate
i=int(l/(2*dx))
j=int(h/(2*dy))
print('\n   temperature at the middle of the plate: {:5.3f}\n'.format(T[i,j])) 

# generating the grid
x=np.linspace(0,l,m+1)   
y=np.linspace(0,h,n+1)   
Y,X=np.meshgrid(y,x)

# temperature at the vertical middle line
fig1=plt.figure(2)
plt.plot(T[i,:],y)
plt.xlabel('T($^oC$)');plt.ylabel('y(m)')
plt.ylim(0,h)

# temperature contour
fig2=plt.figure(3)
plt.contourf(X,Y,T,50,cmap='jet')
plt.axes().set_aspect('equal')
plt.xlabel('x(m)');plt.ylabel('y(m)')
 
# temperature surface   
fig3=plt.figure(4)
ax=plt.axes(projection='3d')
ax.plot_surface(X,Y,T,cmap='jet')
ax.set_xlabel('x(m)');ax.set_ylabel('y(m)');ax.set_zlabel('T($^oC$)')
plt.show()