from matplotlib.pyplot import *
from numpy import *
import sys

spins = sys.argv[1]
temps = [1, 2.4]
#temp = sys.argv[2] 


infile1 = open("Output/Output_"+spins+"_1.txt")
infile2 = open("Output/Output_"+spins+"_2.4.txt")

infile1.readline()
infile2.readline()

t1 = []
E1 = []
E1mean = []
M1 = []
M1mean = []
a1 = []

t2 = []
E2 = []
E2mean = []
M2 = []
M2mean = []
a2 = []

for line in infile1: 
	words = line.split()
	t1.append(float(words[0]))
	E1.append(float(words[1]))
	E1mean.append(float(words[2]))
	M1.append(float(words[3]))
	M1mean.append(float(words[4]))
	a1.append(float(words[5]))

for line in infile2: 
	words = line.split()
	t2.append(float(words[0]))
	E2.append(float(words[1]))
	E2mean.append(float(words[2]))
	M2.append(float(words[3]))
	M2mean.append(float(words[4]))
	a2.append(float(words[5]))

infile1.close()
infile2.close()

t1 = array(t1) 
E1 = array(E1) 
E1mean = array(E1mean)
M1 = array(M1) 
M1mean = array(M1mean)
a1 = array(a1)

t2 = array(t2) 
E2 = array(E2) 
E2mean = array(E2mean)
M2 = array(M2) 
M2mean = array(M2mean)
a2 = array(a2)

"""
figure() 
plot(t1, E1mean) 
xlabel('# MC cycles')
ylabel('Energy')
title('L=20, T=1')
#xlim(0,50000)
#ylim(-525, -475)
#show() 
savefig("Output/Energy_"+spins+"_1.png")

figure() 
plot(t1, M1mean) 
xlabel('# MC cycles')
ylabel('Magnetization')
title('L=20, T=1')
#xlim(0,10000)
#ylim(399,400)
#show() 
savefig("Output/Magnetization_"+spins+"_1.png")

figure() 
plot(t2, E2mean) 
xlabel('# MC cycles')
ylabel('Energy')
title('L=20, T=2.4')
#xlim(0,50000)
#ylim(-525, -475)
#show() 
savefig("Output/Energy_"+spins+"_2.4.png")

figure() 
plot(t2, M2mean) 
xlabel('# MC cycles')
ylabel('Magnetization')
title('L=20, T=2.4')
#xlim(0,10000)
#ylim(399,400)
#show() 
savefig("Output/Magnetization_"+spins+"_2.4.png")

figure()
plot(t1,a1) 
plot(t1,a2)
legend(['T = 1', 'T=2.4'])
xlabel('# MC cycles')
ylabel('Accepted configurations')
#show()
savefig("Output/Accepted_configurations_"+spins+".png")

data = E1[30000:-1]
binwidth = 0.001

figure()
hist(data,bins = arange(min(data), max(data)+binwidth, binwidth))
savefig("Output/Hist_energies.png")
"""

Energies1 = unique(E1[200000:-1])
Energies2 = unique(E2[200000:-1])
E1 = E1[200000:-1]
E2 = E2[200000:-1]
t = t1[200000:-1]

prob_dist1 = {}
prob_dist2 = {}

for i in range(len(Energies1)): 
	prob_dist1[str(Energies1[i])] = 0 

for i in range(len(Energies2)): 
	prob_dist2[str(Energies2[i])] = 0 

norm = len(t)
for i in range(len(t)): 
	prob_dist1[str(E1[i])] += 1.0/norm
	prob_dist2[str(E2[i])] += 1.0/norm

p1 = []
#print 'PDF for T = 1.0:'
for E in Energies1:
	p = prob_dist1[str(E)]
	p1.append(p) 
	#print 'P(%d) = %f' %(E,p)

p2 = []
#print 'PDF for T = 2.4'
for E in Energies2:
	p = prob_dist2[str(E)] 
	p2.append(p)
	#print 'P(%d) = %f' %(E,p)

p1 = array(p1) 
p2 = array(p2)

figure()
plot(Energies1,p1)
xlabel('E')
ylabel('P(E)')
savefig('Output/Energy_PDF_T=1.png')

figure()
plot(Energies2,p2)
xlabel('E')
ylabel('P(E)')
savefig('Output/Energy_PDF_T=2.4.png')






