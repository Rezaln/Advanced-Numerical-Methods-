"""
This code was my university project 
Code by: Reza Lotfi Navaei

"""

#===================================
#           Gauss_Seidel
#===================================
import math
import numpy as np
import matplotlib.pyplot as plt
import time


#============================= main program ============================
Gs_start = time.time()

n=60                          # number of divisions
pi=math.acos(-1)              # pi number
l=2*pi                        # length of the domain
dx=l/n
w=1.                          # over relaxation factor
cr=1e-10                      # convergence criteria

# matrix definition
x=np.zeros(n+1);y=np.zeros(n+1);y_old=np.zeros(n+1)

# boundary conditions
y[0]=7.0
y[n]=0

# generating the grid
for i in range(n+1):
    x[i]=i*dx
    
    
#=============================== solution ==============================
k=0
err=1
while err>cr:
    k=k+1
    
    err=0.
    y_old[:]=y[:]
    for i in range(1,n):
        y[i]=(1-w)*y[i]+(w/2)*(y[i-1]+y[i+1])
#        y[i]=(1-w)*y[i]+(w/(2-3*dx**2))*(y[i-1]+y[i+1])
        err=err+abs(y[i]-y_old[i])/(n-1)
       
    # printing iteration and error
    print('{:6d}{:15.2e}'.format(k,err))
    


#=============================== results ===============================
Gs_end = time.time()
Gs_elapsed_time = Gs_end - Gs_start

print ('Gaus-seidel elapsed time : ', Gs_elapsed_time)

# ploting the result
plt.figure()
plt.plot(x,y)
plt.xlabel('x');plt.ylabel('y')
plt.xlim(0,2*pi)