from numpy import array, zeros 
from matplotlib.pyplot import *
import sys 
from math import sqrt,pi

font = {'size':16}
matplotlib.rc('font', **font)

b = float(sys.argv[1]) 
#v=1

infile = open("Output/Binary_VV_beta=%g.txt" %b, 'r')
infile.readline() 

x = []
y = []

for line in infile: 
	words = line.split()
	x.append(float(words[3]))
	y.append(float(words[4]))

infile.close() 

x = array(x)
y = array(y)


figure() 
plot(x,y)
title(r'$\beta=%g$'%b) 
xlabel(r'$x (AU)$')
ylabel(r'$y (AU)$')
axis('equal')
xlim(min(x)-0.5,max(x)+0.5)
ylim(min(y)-0.5, max(y)+0.5)
#axis('equal')
savefig('Output/Modified_gravity_beta=%g.png' %b) 
#show()

