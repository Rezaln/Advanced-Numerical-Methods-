"""
This code was my university project 
Code by: Reza Lotfi Navaei

"""

#===================================
#               TDMA
#===================================
import math
import numpy as np
import matplotlib.pyplot as plt
import time

#========================= function definition =========================
TDMA_start = time.time()

def TDMA(a,b,c,r,s):

    n=len(a)
    # forward substitution
    c[0]=c[0]/b[0]
    r[0]=r[0]/b[0]
    for i in range(1,n):
        c[i]=c[i]/(b[i]-a[i]*c[i-1])
        r[i]=(r[i]-a[i]*r[i-1])/(b[i]-a[i]*c[i-1])
    
    # backward substitution
    s[n-1]=r[n-1]
    for i in reversed(range(n-1)):
        s[i]=r[i]-s[i+1]*c[i]
    


#============================= main program ============================
n=50                         # number of divisions
pi=math.acos(-1)              # pi number
l=2*pi                        # length of the domain
dx=l/n

# matrix definition
x=np.zeros(n+1);y=np.zeros(n+1);y_exact=np.zeros(n+1)
a=np.zeros(n-1);b=np.zeros(n-1)
c=np.zeros(n-1);r=np.zeros(n-1);s=np.zeros(n-1)

# boundary conditions
y[0]=7.
y[n]=0.

# generating the grid
for i in range(n+1):
    x[i]=i*dx
    

#=============================== solution ==============================
for i in range(1,n):
    a[i-1]=1.
    b[i-1]=-2+3*dx**2
    c[i-1]=1.
    r[i-1]=0.
    
i=1
r[i-1]=r[i-1]-a[i-1]*y[0]
i=n-1
r[i-1]=r[i-1]-c[i-1]*y[n]

# solving the system of equations
TDMA(a,b,c,r,s)

for i in range(1,n):
    y[i]=s[i-1]


#=============================== results ===============================
#Time
TDMA_end = time.time()
TDMA_elapsed_time = TDMA_end - TDMA_start

print ( "TDMA elapsed time : " , TDMA_elapsed_time)

# calculating the exact solution
for i in range(n+1):
    y_exact[i]=7./math.sin(-2*pi*math.sqrt(3))*math.sin( math.sqrt(3)*x[i]-2*pi*math.sqrt(3) )
    
# ploting the numerical solution
plt.figure()
plt.plot(x,y,'b.')
plt.xlabel('x');plt.ylabel('y')
plt.xlim(0,2*pi)
plt.title('TDMA solution')

# ploting the exact solution
plt.figure()
plt.plot(x,y_exact)
plt.xlabel('x');plt.ylabel('y')
plt.xlim(0,2*pi)
plt.title('exact solution')

# ploting the error
plt.figure()
plt.plot(x,y-y_exact,'r')
plt.xlabel('x');plt.ylabel('y')
plt.xlim(0,2*pi)
plt.title('error')