"""
This code was my university project 
Code by: Reza Lotfi Navaei

"""

#===================================
#          Modified_Euler
#===================================
import math
import numpy as np
import matplotlib.pyplot as plt
import time

#========================= function definition =========================
ME_start = time.time()

def f(x,y,z):
    return -3*y

def g(x,y,z):
    return z    


#============================= main program ============================
n=100                         # number of divisions
pi=math.acos(-1)              # pi number
l=2*pi                        # length of the domain
dx=l/n
inc=1e-4                      # shooting method's increment
cr=1e-5                       # convergence criteria

# matrix definition
x=np.zeros(n+1);y=np.zeros(n+1);z=np.zeros(n+1)

# boundary conditions
y[0]=7.
z[0]=0.                       # initial guess

# generating the grid
for i in range(n+1):
    x[i]=i*dx
    
    
#=============================== solution ==============================
k=0
err=1.
while err>cr: 
    
    k=k+1
    z[0]=z[0]-inc

    for i in range(n):
        
        k1=f( x[i]    ,y[i]        ,z[i])
        m1=g( x[i]    ,y[i]        ,z[i])
        
        k2=f( x[i+1]  ,y[i]+dx*m1  ,z[i]+dx*k1)
        m2=g( x[i+1]  ,y[i]+dx*m1  ,z[i]+dx*k1)
        
        
        z[i+1]=z[i]+(dx/2.)*( k1+k2 )
        y[i+1]=y[i]+(dx/2.)*( m1+m2 )
        
    err=-y[n]-0
    
    # printing iteration, error, first derivatve at the left boundary and  calculated right boundary value
    print('{:6d}{:15.2e}{:10.5f}{:10.5f}'.format(k,err,z[0],y[n]))



    
#=============================== result ================================
ME_end = time.time()
ME_elapsed = ME_end - ME_start

print ("Modified euler elapsed time : ",ME_elapsed)

# ploting the result
plt.figure()
plt.plot(x,y)
plt.xlabel('x');plt.ylabel('y')
plt.xlim(0,l)