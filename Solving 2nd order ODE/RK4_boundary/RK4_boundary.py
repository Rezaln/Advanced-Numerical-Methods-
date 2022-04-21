"""
This code was my university project 
Code by: Reza Lotfi Navaei

"""

#===================================
#               RK4
#===================================
import math
import numpy as np
import matplotlib.pyplot as plt
import time

#========================= function definition =========================
RK4_start = time.time()

def f(x,y,z):
    return -3*y

def g(x,y,z):
    return z    



#============================= main program ============================
n=50                         # number of divisions
pi=math.acos(-1)              # pi number
l=2*pi                        # length of the domain
dx=l/n

# matrix definition
x=np.zeros(n+1);y=np.zeros(n+1)
z=np.zeros(n+1);y_exact=np.zeros(n+1);y_f=np.zeros((2,n+1))

# boundary conditions
y[0]=7.

# generating the grid
for i in range(n+1):
    x[i]=i*dx
    
    
#=============================== solution ==============================
k=-1
for z[0] in [-1,-2]:          # initial guess
    k=k+1

    for i in range(n):
        
        k1=f( x[i]         ,y[i]            ,z[i])
        m1=g( x[i]         ,y[i]            ,z[i])
        
        k2=f( x[i]+0.5*dx  ,y[i]+0.5*dx*m1  ,z[i]+0.5*dx*k1)
        m2=g( x[i]+0.5*dx  ,y[i]+0.5*dx*m1  ,z[i]+0.5*dx*k1)
        
        k3=f( x[i]+0.5*dx  ,y[i]+0.5*dx*m2  ,z[i]+0.5*dx*k2)
        m3=g( x[i]+0.5*dx  ,y[i]+0.5*dx*m2  ,z[i]+0.5*dx*k2)
        
        k4=f( x[i+1]       ,y[i]+dx*m3      ,z[i]+dx*k3)
        m4=g( x[i+1]       ,y[i]+dx*m3      ,z[i]+dx*k3)
        
        
        z[i+1]=z[i]+(dx/6.)*( k1+2*(k2+k3)+k4 )
        y[i+1]=y[i]+(dx/6.)*( m1+2*(m2+m3)+m4 )
        
    # saving solutions for interpolating later
    for i in range(n+1):
        y_f[k,i]=y[i]


# coefficients of interpolation        
c1=y_f[1,n]/(y_f[1,n]-y_f[0,n])   
c2=y_f[0,n]/(y_f[0,n]-y_f[1,n])   

# interpolating between solutions
for i in range(i+1):
    y_f[0,i]=c1*y_f[0,i]+c2*y_f[1,i]
    
#=============================== Time ==================================
RK4_end  = time.time()

RK4_elapsed = RK4_end - RK4_start

print ( 'RK4 elapsed : ' , RK4_elapsed )

#=============================== result ================================
# calculating the exact solution
for i in range(n+1):
    y_exact[i]=7./math.sin(-2*pi*math.sqrt(3))*math.sin( math.sqrt(3)*x[i]-2*pi*math.sqrt(3) )
    
# ploting the numerical solution
plt.figure()
plt.plot(x,y_f[0,:],'b.')
plt.xlabel('x');plt.ylabel('y')
plt.xlim(0,2*pi)
plt.title('RK4 solution')

# ploting the exact solution
plt.figure()
plt.plot(x,y_exact)
plt.xlabel('x');plt.ylabel('y')
plt.xlim(0,2*pi)
plt.title('exact solution')

# ploting the error
plt.figure()
plt.plot(x,y_f[0,:]-y_exact,'r')
plt.xlabel('x');plt.ylabel('y')
plt.xlim(0,2*pi)
plt.title('error')