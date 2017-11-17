#from math import *
from numpy import array, zeros 
from matplotlib.pyplot import *
import sys 
from math import sqrt

font = {'size':16}
matplotlib.rc('font', **font)

problem = sys.argv[1]  
mJ = int(sys.argv[2]) 

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


x_max = max(x['earth'])
y_max = max(y['earth'])
x_min = min(x['earth'])
y_min = min(y['earth'])

print "x max = %f" %x_max 
print "y max = %f" %y_max
print "x min = %f" %x_min 
print "y min = %f" %y_min

x_s = array([0.0])
y_s = array([0.0])
z_s = array([0.0]) 

n = len(x['earth'])

if problem == "SolarSystem": 
	from mpl_toolkits.mplot3d import Axes3D 
	from matplotlib.font_manager import FontProperties
	fontP = FontProperties() 
	fontP.set_size('small')
	figure()
	ax = axes(projection = '3d') 
	for planet in planets: 
		ax.plot(x[planet], y[planet], z[planet], label = planet)
	#plot(x_s,y_s,'ro')
	ax.legend(prop=fontP)
	ax.set_xlabel(r'$x (AU)$')
	ax.set_ylabel(r'$y (AU)$')
	#ax.set_zlabel(r'$z (AU)$')
	savefig('Output/SolarSystem_3D.png')
	#show()

figure()
for planet in planets: 
	#if planet == "sun": 
	#	continue
	plot(x[planet], y[planet], label = planet)
legend()
xlabel(r'$x (AU)$')
ylabel(r'$y (AU)$')
#if problem == "Binary": 
#axis([-1.5, 1.5, -1.5, 1.5])
axis('equal')
savefig('Output/'+problem+'_n=%d.png' %n)
#show()
 
	
