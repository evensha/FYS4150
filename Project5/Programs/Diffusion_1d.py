from matplotlib.pyplot import *
from numpy import *
import sys

infile_FE = open('Output/Output_FE.txt', 'r')
infile_BE = open('Output/Output_BE.txt', 'r')
infile_CN = open('Output/Output_CN.txt', 'r')

u1_FE = []
u2_FE = []

u1_BE = []
u2_BE = []

u1_CN = []
u2_CN = []

for line in infile_FE: 
	words = line.split()
	u1_FE.append(float(words[0]))
	u2_FE.append(float(words[1])) 

for line in infile_BE: 
	words = line.split()
	u1_BE.append(float(words[0]))
	u2_BE.append(float(words[1])) 

for line in infile_CN: 
	words = line.split() 
	u1_CN.append(float(words[0])) 
	u2_CN.append(float(words[1]))

N = len(u1_FE)

#print N

x = linspace(0,1,N)
u_ana_1 = zeros(N)
u_ana_2 = zeros(N) 
u_ana_1[N-1] = 1.0
u_ana_2[N-1] = 1.0

t1 = 0.5*(1/float(N))**2*10
t2 = 0.5*(1/float(N))**2*100
print t1 
print t2


for i in range(1,len(x)-1):
	for n in range(1,25): 
		u_ana_1[i] += (-1)**(n-1)/float(n)*sin(n*pi*x[i])*exp(-(n*pi)**2*t1)			
		u_ana_2[i] += (-1)**(n-1)/float(n)*sin(n*pi*x[i])*exp(-(n*pi)**2*t2)	
	u_ana_1[i] = x[i] - 2/pi*u_ana_1[i]
	u_ana_2[i] = x[i] - 2/pi*u_ana_2[i]

print u_ana_1
#print u_ana_2
"""
figure(); 
plot(x,u_ana_1)
plot(x,u_ana_2)
show()
"""
u1_FE = array(u1_FE)
u2_FE = array(u2_FE)
u1_BE = array(u1_BE)
u2_BE = array(u2_BE)
u1_CN = array(u1_CN)
u2_CN = array(u2_CN) 

#print len(u_ana_1) 
#print len(u1_CN)
#print x

figure(); 
plot(x, u1_FE) 
plot(x, u1_BE) 
plot(x, u1_CN)
plot(x, u_ana_1)
legend(['FE', 'BE', 'CN', 'Ana'], loc = 2)
xlabel('x')
ylabel('u(x)')
title(r'$\Delta x = 1/10, t=%.2f$' %t1 )
savefig('Output/Diffusion_1d_t1.png')
show()

figure(); 
plot(x, u2_FE)
plot(x, u2_BE)
plot(x, u2_CN)
plot(x, u_ana_2)
legend(['FE', 'BE', 'CN', 'Ana'])
xlabel('x')
ylabel('u(x)')
title(r'$\Delta x = 1/10, t=%.2f$' %t2 )
savefig('Output/Diffusion_1d_t2.png')
#show()


