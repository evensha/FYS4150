from matplotlib.pyplot import *
from numpy import *
import sys

font = {'size':18}
matplotlib.rc('font', **font)

N =  int(sys.argv[1])

infile_FE = open('Output/Output_FE_%d.txt' %N, 'r')
infile_BE = open('Output/Output_BE_%d.txt' %N, 'r')
infile_CN = open('Output/Output_CN_%d.txt' %N, 'r')

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


x = linspace(0,1,N+1)
u_ana_1 = zeros(N+1)
u_ana_2 = zeros(N+1) 
u_ana_1[N] = 1.0
u_ana_2[N] = 1.0

t1 = 0.05
t2 = 0.5
#print t1 
#print t2


for i in range(1,len(x)-1):
	for n in range(1,25): 
		u_ana_1[i] += (-1)**(n-1)/float(n)*sin(n*pi*x[i])*exp(-(n*pi)**2*t1)			
		u_ana_2[i] += (-1)**(n-1)/float(n)*sin(n*pi*x[i])*exp(-(n*pi)**2*t2)	
	u_ana_1[i] = x[i] - 2/pi*u_ana_1[i]
	u_ana_2[i] = x[i] - 2/pi*u_ana_2[i]

#print u_ana_1
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
xlabel(r'$x$')
ylabel(r'$u(x)$')
title(r'$\Delta x = 1/%d, t=%.2f$' %(N+1,t1) )
savefig('Output/Diffusion_1d_%d_t1.png' %N, bbox_inches = 'tight')
#show()

figure(); 
plot(x, u2_FE)
plot(x, u2_BE)
plot(x, u2_CN)
plot(x, u_ana_2)
legend(['FE', 'BE', 'CN', 'Ana'], loc = 2)
xlabel(r'$x$')
ylabel(r'$u(x)$')
title(r'$\Delta x = 1/%d, t=%.2f$' %(N+1,t2) )
savefig('Output/Diffusion_1d_%d_t2.png' %N, bbox_inches = 'tight')
#show()
"""
figure(); 
plot(x, abs(u1_FE-u_ana_1))
plot(x, abs(u1_BE-u_ana_1))
plot(x, abs(u1_CN-u_ana_1))
xlabel(r'$x$')
ylabel(r'$\Delta u$')
legend(['FE', 'BE', 'CN'], loc = 2)
savefig('Output/Delta_t1_1d_%d.png' %N, bbox_inches = 'tight')
#show()

figure(); 
plot(x, abs(u2_FE-u_ana_2))
plot(x, abs(u2_BE-u_ana_2))
plot(x, abs(u2_CN-u_ana_2))
xlabel(r'$x$')
ylabel(r'$\Delta u$')
legend(['FE', 'BE', 'CN'], loc = 2)
savefig('Output/Delta_t2_1d_%d.png' %N, bbox_inches = 'tight')
#show()
"""
