from math import * 
from numpy import * 
from matplotlib.pyplot import *
import sys
import time 

alg = sys.argv[1];  # indicates if we want to run the general (g) or the simplified (s) algorithm
n = int( sys.argv[2] ); # matrix dimension 

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
t0 = time.clock() 


# Forward substitution:  

if alg == 'g':  
	for i in range(2,n+1):  # start loop on b[2] and f[2]  
		#print i
		b[i] = b[i] - a[i-1]*c[i-1]/b[i-1]   
		f[i] = f[i] - a[i-1]*f[i-1]/b[i-1]
		fl_ops += 6


# Simplified forward substitution: 

if alg == 's': 
	for i in range(2,n+1): 
		b[i] = (i+1)/float(i)    
		f[i] = f[i] + (i-1)/float(i)*f[i-1] 
		fl_ops += 4 


# Backward substitution:  

v = zeros(n+2)   
v[0]= 0; v[n+1]= 0 # boundary values 
v[n] = f[n]/b[n] # n'th element of u  

if alg == 'g':
	for i in range(n,1,-1):  # loop backwards from the second last element of u    
		#print i
		v[i-1] = (f[i-1] - c[i-1]*v[i])/b[i-1]  
		fl_ops += 3


# Simplified backward substitution:

if alg == 's': 
	for i in range(n,1,-1): 
		#print i
		v[i-1] = (i-1)/float(i)*(f[i-1]+v[i]) 
		fl_ops += 3


# Floating point operations and CPU time: 

t1 = time.clock() 
cpu_time = t1- t0 
print 'CPU time: %f' %cpu_time

print 'Floting point operations: %d' %fl_ops


# Analytical solution 

u = [1-(1-exp(-10))*x[i] - exp(-10*x[i]) for i in range(len(x))]


# Plot numerical and analytical solution together

plot(x,v,'r')
plot(x,u,'b')
legend(['Numerical','Analytical'])
xlabel('x')
ylabel('u(x)')
show()


# Errors 

error = zeros(n+2) 

for i in range(1,n+1): 
	error[i] = abs((v[i]-u[i])/u[i])

max_error = max(error)
log_error = log10(max_error) 

print 'Error with n=%d: %f' %(n,log_error)

"""
A = array(zeros((n,n))) 

for i in range(n): 
	for j in range(n): 
		if i == j: 
			A[i][j] = 2 
		if i == j-1 or j == i-1: 
			A[i][j] = -1 

f_1 = [f[i] for i in range(1,n+1)]
print A
from numpy.linalg import solve 

#print f 
#print f_1

u_1 = solve(A,f_1)
print u_1
print u_ana
"""


