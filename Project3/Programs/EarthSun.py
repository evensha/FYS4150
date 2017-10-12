#from math import *
from numpy import array 
from matplotlib.pyplot import * 

infile1 = open("Output/EarthSun_VV.txt", 'r')
infile2 = open("Output/EarthSun_FE.txt", 'r')

x1 = []
y1 = []
x2 = []
y2 = []

for line in infile1: 
	words = line.split()
	x1.append(float(words[0]))
	y1.append(float(words[1]))

for line in infile2:
	words = line.split()
	x2.append(float(words[0]))
	y2.append(float(words[1]))


infile1.close()
infile2.close()

#print x
#print y

x1 = array(x1)
y1 = array(y1)
x2 = array(x2)
y2 = array(y2)

x_s = array([0.0])
y_s = array([0.0]) 

figure() 
plot(x1,y1,'g', label='Verlet')
#plot(x2,y2,'b', label='Euler')
plot(x_s,y_s,'ro')
legend()
xlabel(r'$x$')
ylabel(r'$y$')
#axis('equal')
show()
savefig('Output/EarthSun.png')

 
	
