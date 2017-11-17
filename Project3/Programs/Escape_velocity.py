from numpy import array, zeros 
from matplotlib.pyplot import *
import sys 
from math import sqrt,pi

font = {'size':16}
matplotlib.rc('font', **font)

v = float(sys.argv[1]) 
#v=1

infile = open("Output/Binary_VV.txt", 'r')
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

v1 = 2*v

figure() 
plot(x,y)
title(r'$v=%.1f\pi$ AU/year'%v1) 
xlabel(r'$x (AU)$')
ylabel(r'$y (AU)$')
axis('equal')
xlim(min(x)-0.5,max(x)+0.5)
ylim(min(y)-0.5, max(y)+0.5)
#axis('equal')
savefig('Output/Escape_v=%.1f.png' %v) 
show()


