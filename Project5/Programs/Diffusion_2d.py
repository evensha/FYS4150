from matplotlib.pyplot import *
from numpy import *
import sys

from mpl_toolkits.mplot3d import Axes3D

infile = open('Output/Diffusion_2d.txt', 'r')

n = 10 

#x,y = meshgrid(x,y)

#z =  sqrt(x +y)

u = zeros((n,n))
i = 0

#x = arange(0,1.01,1/9.0)
#y = arange(0,1.01,1/9.0) 
x = linspace(0,1,n)
y = linspace(0,1,n) 

#x,y = meshgrid(x,y)

for line in infile:
	words = line.split() 
	for j in range(0,n):  
		u[i,j] = float(words[j])
	i+= 1

#print u
#z = sqrt(x+y)

fig = figure()
ax = fig.add_subplot(111, projection = '3d')

surf = ax.plot_surface(x,y,u)
show()



