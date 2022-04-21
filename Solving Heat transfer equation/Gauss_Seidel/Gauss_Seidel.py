#===================================
#           Gauss_Seidel
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
alpha=23.1*1e-6               # diffusivity of the iron
S_T=5e-3                      # heat source
w=1                         # over relaxation factor
cr=1e-6                       # convergence criteria


# initial guess
T=np.zeros((m+1,n+1));TT=np.zeros((m+1,n+1));

# boundary conditions
for i in range(m+1):
    T[i,0]=20
    T[i,n]=40

for j in range(n+1):
    T[m,j]=10



#=============================== solution ==============================
k=0
err=1
while err>cr:
    k=k+1
    
    # updating boundary conditions
    for j in range(n+1):
        T[0,j]=T[1,j]
        
    TT[:,:]=T[:,:]
            
    # Gauss_Seidel iteration
    err=0
    for i in range(1,m):
        for j in range(1,n):
            T[i,j]=(1-w)*T[i,j]+(w/(2*(1+be**2)))*( T[i,j-1]+T[i,j+1] +be**2*( T[i-1,j]+T[i+1,j] )  +S_T/alpha*dy**2 )
            err=err+math.sqrt( (T[i,j]-TT[i,j])**2/((m-1)*(n-1)) )
     
        
    plt.scatter(k,err)    
    plt.xlabel('iteration');plt.ylabel('norm of the error')
    plt.title('convergence history')
    plt.yscale("log")
    plt.pause(0.005)
    
    # printing iteration and second norm of the error
    print('{:6d}{:15.2e}'.format(k,err))
    
    

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