#===================================
#                FTCS
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
dt=5                          # time step
Tf=30000                      # final time
alpha=23.1*1e-6                 # diffusivity of the iron
S_T=5e-3                      # heat source
cr=1e-3                       # convergence criteria
print(alpha)
Fo_x=alpha*dt/dx**2
Fo_y=alpha*dt/dy**2

if Fo_x>0.25 or 0.25<Fo_y:
    raise Exception('time step is large, reduce it!')

# initial condition
T=np.zeros((m+1,n+1));T_old=np.zeros((m+1,n+1));

# boundary conditions
for i in range(m+1):
    T[i,0]=20
    T[i,n]=40

for j in range(n+1):
    T[m,j]=10



#=============================== solution ==============================
t=0
err=1
while err>cr and t<Tf:
    t=t+dt
    
    # updating boundary conditions
    for j in range(n+1):
        T[0,j]=T[1,j]
        
    T_old[:,:]=T[:,:]
        
    # FTCS loop
    err=0
    for i in range(1,m):
        for j in range(1,n):
            T[i,j]=T_old[i,j] +Fo_x*( T_old[i-1,j]-2*T_old[i,j]+T_old[i+1,j] )  +Fo_y*( T_old[i,j-1]-2*T_old[i,j]+T_old[i,j+1] )   +dt*S_T
            err=err+math.sqrt( (T[i,j]-T_old[i,j])**2/((m-1)*(n-1)) )
     
        
#    plt.scatter(t,err)    
#    plt.xlabel('time(s)');plt.ylabel('norm of the error')
#    plt.title('convergence history')
#    plt.yscale("log")
#    plt.pause(0.0001)
    
    # printing time and second norm of the error
    print('{:10.2f}{:15.2e}'.format(t,err))
    
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