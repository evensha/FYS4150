from matplotlib.pyplot import *
from numpy import *
import sys
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter

infile = open('Output/Diffusion_2d.txt', 'r')

n = 100

u = zeros((n,n))
u_ana = zeros((n,n))
i = 0

x = linspace(0,1,n)
y = linspace(0,1,n) 

#print x

for line in infile:
	words = line.split() 
	for j in range(0,n):  
		u[i,j] = float(words[j])
	i+= 1

for i in range(0,n): 
	u_ana[i,n-1] = i/float(n-1)
	u_ana[n-1,i] = i/float(n-1)

#print u_ana

t1 = 0.05

for i in range(1,len(x)-1):
	for j in range(1,len(y)-1): 
		for n in range(1,15):
			for m in range(1,15):
				#u_ana[i,j] += 1.0  
				u_ana[i,j] += (-1)**(n+m-2)/(float(n*m))*sin(n*pi*x[i])*sin(m*pi*y[j])*exp(-(n**2 + m**2)*pi**2*t1)
		u_ana[i,j] = x[i]*y[j] - 4/pi**2*u_ana[i,j]

#print u_ana

for i in range(0,n): 
	for j in range(0,n): 
		if u[i,j]<0: 
			print "Negative element!"

x,y = meshgrid(x,y)


fig = figure()
ax = fig.gca(projection = '3d')

surf = ax.plot_surface(x,y,u, rstride = 5, cstride = 5, cmap=cm.coolwarm,linewidth = 0,antialiased = False)
fig.savefig('Output/Diffusion_2d.png')
#show()

fig1 = figure()
ax1 = fig1.gca(projection = '3d')

surf = ax1.plot_surface(x,y,abs(u-u_ana),rstride = 5, cstride = 5, cmap=cm.coolwarm,linewidth = 0,antialiased = False)
fig1.savefig('Output/Delta_2d.png')
#show()



