from numpy import array, zeros 
from matplotlib.pyplot import *
import sys 
from math import sqrt,pi

#v = 2*pi*float(sys.argv[1]) 
v=1

infile = open("Output/Binary_VV.txt", 'r')
infile.readline() 

x = []
y = []

for line in infile: 
	words = line.split() 
	x.append(float(words[3])
	y.append(float(words[4])

infile.close() 

x = array(x)
y = array(y)

figure() 
plot(x,y) 
xlabel(r'$x$')
ylabel(r'$y$')
axis('equal')
savefig('Output/FE_vs_VV_v=%f.png' %v) 
show()


