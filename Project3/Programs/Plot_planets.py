#from math import *
from numpy import array 
from matplotlib.pyplot import *
import sys 

problem = sys.argv[1]  

planets = []

if problem == "SolarSystem": 
	infile = open("Output/SolarSystem_VV.txt", 'r')
	planets = ['sun','mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptun', 'pluto']
if problem == "ThreeBody": 
	infile = open("Output/ThreeBody_VV.txt", 'r') 
	planets = ['sun','earth', 'jupiter']
if problem == "Binary": 
	infile = open("Output/Binary_VV.txt", 'r')
	planets = ['sun','earth']


infile.readline()

x = {}
y =  {}

x_place =  {}
y_place =  {}

i = 0
for planet in planets: 
	#print planet
	x[planet] = []
	y[planet] = []
	x_place[planet] = i 
	y_place[planet] = i+1
	i+= 2	
	

#print x_place 
#print y_place


#print x.keys()

for line in infile: 
	words = line.split()
	for planet in x.keys(): 
		x[planet].append(float(words[x_place[planet]]))
		y[planet].append(float(words[y_place[planet]]))

infile.close()
	
for planet in planets: 
	x[planet] = array(x[planet])
	y[planet] = array(y[planet])

x_s = array([0.0])
y_s = array([0.0])


figure() 
for planet in planets: 
	plot(x[planet], y[planet], label=planet)
#plot(x_s,y_s,'ro')
legend()
xlabel(r'$x$')
ylabel(r'$y$')
#if problem == "Binary": 
#axis([-1.5, 1.5, -1.5, 1.5])
axis('equal')
show()
savefig('Output/'+problem+'.png')

 
	
