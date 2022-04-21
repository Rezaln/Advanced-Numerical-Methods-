
"""
This code was my university project 
Code by: Reza Lotfi Navaei

"""
#===================================
#               PDMA
#===================================
import math
import numpy as np
import matplotlib.pyplot as plt
import time

#========================= function definition =========================
PDMA_start = time.time()


def PDMA(e,a,b,c,f,r,s):

    n=len(e)
    # first forward substitution
    for i in range(2,n):
        a[i]=a[i]-e[i]*b[i-1]/a[i-1]
        b[i]=b[i]-e[i]*c[i-1]/a[i-1]
        c[i]=c[i]-e[i]*f[i-1]/a[i-1]
        r[i]=r[i]-e[i]*r[i-1]/a[i-1]
    
    # second forward substitution
    c[0]=c[0]/b[0]
    f[0]=f[0]/b[0]
    r[0]=r[0]/b[0]
    for i in range(1,n):
        c[i]=(c[i]-a[i]*f[i-1])/(b[i]-a[i]*c[i-1]);
        f[i]=f[i]/(b[i]-a[i]*c[i-1]);
        r[i]=(r[i]-a[i]*r[i-1])/(b[i]-a[i]*c[i-1]);
        
    # backward substitution
    s[n-1]=r[n-1]
    s[n-2]=r[n-2]-c[n-2]*s[n-1]
    for i in reversed(range(n-2)):
        s[i]=r[i]-s[i+1]*c[i]-s[i+2]*f[i]
    


#============================= main program ============================
n=50                          # number of divisions
pi=math.acos(-1)              # pi number
l=2*pi                        # lengtho f the domain
dx=l/n

# matrix definition
x=np.zeros(n+1);y=np.zeros(n+1);y_exact=np.zeros(n+1)
e=np.zeros(n-1);a=np.zeros(n-1)
b=np.zeros(n-1);c=np.zeros(n-1)
f=np.zeros(n-1);r=np.zeros(n-1);s=np.zeros(n-1)

# boundary conditions
y[0]=7
y[n]=0

# generating the grid
for i in range(n+1):
    x[i]=i*dx
    
    
#=============================== solution ==============================
for i in range(1,n):
    e[i-1]=-1./12
    a[i-1]=16./12
    b[i-1]=-30./12+3*dx**2
    c[i-1]=16./12
    f[i-1]=-1./12
    r[i-1]=0

# imposing boundary conditions for third and second to last nodes
i=2
r[i-1]=r[i-1]-e[i-1]*y[0]
i=n-2
r[i-1]=r[i-1]-f[i-1]*y[n]


# coefficients and right hand side matrix for the second node      
i=1
a[i-1]=1
b[i-1]=-2
c[i-1]=1
f[i-1]=0
r[i-1]=r[i-1]-a[i-1]*y[0]

# coefficients and right hand side matrix for the one to last node     
i=n-1
a[i-1]=1
b[i-1]=-2
c[i-1]=1
e[i-1]=0
r[i-1]=r[i-1]-c[i-1]*y[n]

# solving the system of equations
PDMA(e,a,b,c,f,r,s)

for i in range(1,n):
    y[i]=s[i-1]


#=============================== results ===============================
#Time
PDMA_end = time.time()
PDMA_elapsed_time = PDMA_end - PDMA_start

print (' PDMA elapsed time : ', PDMA_elapsed_time)

# calculating the exact solution
for i in range(n+1):
    y_exact[i]=7./math.sin(-2*pi*math.sqrt(3))*math.sin( math.sqrt(3)*x[i]-2*pi*math.sqrt(3) )
    
# ploting the numerical solution
plt.figure()
plt.plot(x,y,'b.')
plt.xlabel('x');plt.ylabel('y')
plt.xlim(0,2*pi)
plt.title('PDMA solution')

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