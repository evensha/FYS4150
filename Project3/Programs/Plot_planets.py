#from math import *
from numpy import array, zeros 
from matplotlib.pyplot import *
import sys 
from math import sqrt

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
z = {}

x_place =  {}
y_place =  {}
z_place = {}

i = 0
for planet in planets: 
	#print planet
	x[planet] = []
	y[planet] = []
	z[planet] = []
	x_place[planet] = i 
	y_place[planet] = i+1
	z_place[planet] = i+2
	i+= 3	
	

#print x_place 
#print y_place


#print x.keys()

for line in infile: 
	words = line.split()
	for planet in x.keys(): 
		x[planet].append(float(words[x_place[planet]]))
		y[planet].append(float(words[y_place[planet]]))
		z[planet].append(float(words[z_place[planet]]))

infile.close()
	
for planet in planets: 
	x[planet] = array(x[planet])
	y[planet] = array(y[planet])
	z[planet] = array(z[planet])

r = zeros(len(x['earth']))
#print x['earth'][0]
for i in range(len(x['earth'])): 
	r[i] = sqrt( x['earth'][i]**2 + y['earth'][i]**2 ) 

for i in range(len(r)): 
	if( r[i] > 0.3075 and r[i] < 0.3075001): 
		print r[i]
		print x['earth'][i]
		print y['earth'][i]
		print i
		print '---------'


x_s = array([0.0])
y_s = array([0.0])
z_s = array([0.0]) 
"""
from mpl_toolkits.mplot3d import Axes3D 

figure()
ax = axes(projection = '3d') 
for planet in planets: 
	ax.plot(x[planet], y[planet], z[planet])
#plot(x_s,y_s,'ro')
show()
"""
figure()
for planet in planets: 
	plot(x[planet], y[planet], label = planet)
legend()
xlabel(r'$x$')
ylabel(r'$y$')
#if problem == "Binary": 
#axis([-1.5, 1.5, -1.5, 1.5])
axis('equal')
#show()
savefig('Output/'+problem+'.png')

 
	
