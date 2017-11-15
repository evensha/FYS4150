from matplotlib.pyplot import *
from numpy import *
import sys

spins = sys.argv[1]
temp = sys.argv[2] 


infile = open("Output/Output_"+spins+"_"+temp+".txt")

infile.readline()

t = []
E = []
M = []

for line in infile: 
	words = line.split()
	t.append(float(words[0]))
	E.append(float(words[1]))
	M.append(float(words[2]))

infile.close()

t = array(t) 
E = array(E) 
M = array(M) 

figure() 
plot(t, E) 
xlabel('# MC cycles')
ylabel('Energy')
title('L=20, T=%s' %temp)
#xlim(0,10000)
#ylim(-800, -796)
show() 
savefig("Output/Energy_"+spins+"_"+temp+".png")

figure() 
plot(t, M) 
xlabel('# MC cycles')
ylabel('Magnetization')
title('L=20, T=%s' %temp)
#xlim(0,10000)
#ylim(399,400)
#show() 
savefig("Output/Magnetization_"+spins+"_"+temp+".png")

