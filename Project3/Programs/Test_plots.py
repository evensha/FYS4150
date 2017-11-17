from numpy import array, zeros 
from matplotlib.pyplot import *
import sys 
from math import sqrt

infile1 = open("Output/Binary_FE.txt", 'r')
infile2 = open("Output/Binary_VV.txt", 'r')

infile1.readline() 
infile2.readline() 

x_FE = []; y_FE = []; x_VV = []; y_VV = []; 

for line in infile1: 
	words = line.split() 
	x_FE.append(float(words[3]))
	y_FE.append(float(words[4])) 

for line in infile2: 
	words = line.split() 
	x_VV.append(float(words[3]))
	y_VV.append(float(words[4])) 

infile1.close() 
infile2.close()

x_FE = array(x_FE)
y_FE = array(y_FE)
x_VV = array(x_VV)
y_VV = array(y_VV)

n = len(x_FE)-1 

figure() 
plot(x_FE, y_FE, 'b',  label = 'FE')
plot(x_VV, y_VV, 'r', label = 'VV')
legend()
xlabel(r'$x (AU)$')
ylabel(r'$y (AU)$')
axis('equal')
savefig('Output/FE_vs_VV_n=%d.png' %n) 
show()

