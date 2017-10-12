#from math import *
from numpy import array 
from matplotlib.pyplot import * 

infile1 = open("Output/Planets_VV_3.txt", 'r')
infile2 = open("Output/Planets_VV_1.txt", 'r') 

infile1.readline()
infile2.readline()

x = {'earth':[], 'mars':[], 'jupiter':[]}
y =  {'earth':[], 'mars':[], 'jupiter':[]}

x_place =  {'earth':2, 'mars':4, 'jupiter':6}
y_place =  {'earth':3, 'mars':5, 'jupiter':7}

#print x.keys()

for line in infile1: 
	words = line.split()
	for planet in x.keys(): 
		x[planet].append(float(words[x_place[planet]]))
		y[planet].append(float(words[y_place[planet]]))

	
#print x1[0]

x3 = []
y3 = []

for line in infile2:
	words = line.split()
	x3.append(float(words[2]))
	y3.append(float(words[3]))


infile1.close()
infile2.close()

#print x
#print y

for planet in x.keys(): 
	x[planet] = array(x[planet])
	y[planet] = array(y[planet])

x_s = array([0.0])
y_s = array([0.0]) 

figure() 
plot(x['earth'],y['earth'],'g', label='Earth')
plot(x['mars'],y['mars'],'b', label='Mars')
plot(x['jupiter'], y['jupiter'], 'm', label = 'Jupiter')
#plot(x3,y3, 'm', label = 'Earth (binary)')
plot(x_s,y_s,'ro')
legend()
xlabel(r'$x$')
ylabel(r'$y$')
#axis('equal')
show()
savefig('Output/EarthJupiter.png')

 
	
