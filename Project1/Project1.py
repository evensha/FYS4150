from math import * 
from numpy import * 
from matplotlib.pyplot import *
import sys


n = int( sys.argv[1] ); # matrix dimension 

h = 1.0/(n+1) # step size 


# Define and fill necessary arrays:  

a = zeros(n+2)  # NB: I let all arrays have n+2 elements, i.e. i=0,1,...,n,n+1 
b = zeros(n+2) 
c = zeros(n+2)

for i in range(1,n+1): 
	a[i] = -1.0 
	b[i] = 2.0 
	c[i] = -1.0


x = [i*h for i in range(0,n+2)] 

f = [ h**2*100*exp(-10*x[i]) for i in range(len(x)) ]  


fl_ops = 0; # floting point operations

# Forward substitution:  

for i in range(2,n+1):  # start loop on b[2] and f[2]  
	#print i
	b[i] = b[i] - a[i-1]*c[i-1]/b[i-1]   
	f[i] = f[i] - a[i-1]*f[i-1]/b[i-1]
	fl_ops += 6


# Backward substitution:  

u = zeros(n+2)   
u[0]= 0; u[n+1]= 0 # boundary values 
u[n] = f[n]/b[n] # n'th element of u  

for i in range(n,1,-1):  # loop backwards from the second last element of u    
	#print i
	u[i-1] = (f[i-1] - c[i-1]*u[i])/b[i-1]  
	fl_ops += 3


print 'Floting point operations: %d' %fl_ops

# Analytical solution 

u_ana = [1-(1-exp(-10))*x[i] - exp(-10*x[i]) for i in range(len(x))]


# Plot numerical and analytical solution together

plot(x,u,'r')
plot(x,u_ana,'b')
legend(['Numerical','Analytical'])
xlabel('x')
ylabel('u(x)')
show()

