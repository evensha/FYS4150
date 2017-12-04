from matplotlib.pyplot import *
from numpy import *
import sys
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter

font = {'size':16}
matplotlib.rc('font', **font)

infile1 = open('Output/Diffusion_2d_t1_us.txt', 'r')
infile2 = open('Output/Diffusion_2d_t2_us.txt', 'r')

N = 99

u1 = zeros((N+1,N+1))
u2 = zeros((N+1,N+1))
u_ana_1 = zeros((N+1,N+1))
u_ana_2 = zeros((N+1,N+1))

x = linspace(0,1,N+1)
y = linspace(0,1,N+1) 

#print x
i = 0
for line in infile1:
	words = line.split() 
	for j in range(0,N+1):  
		u1[i,j] = float(words[j])
	i+= 1

i = 0
for line in infile2:
	words = line.split() 
	for j in range(0,N+1):  
		u2[i,j] = float(words[j])
	i+= 1

for i in range(0,N+1): 
	u_ana_1[i,N] = i/float(N)
	u_ana_1[N,i] = i/float(N)
	u_ana_2[i,N] = i/float(N)
	u_ana_2[N,i] = i/float(N)


#print u_ana

t1 = 0.01
t2 = 1.0

for i in range(1,len(x)-1):
	for j in range(1,len(y)-1): 
		for n in range(1,20):
			for m in range(1,20):
				#u_ana[i,j] += 1.0  
				u_ana_1[i,j] += (-1)**(n+m-2)/(float(n*m))*sin(n*pi*x[i])*sin(m*pi*y[j])*exp(-(n**2 + m**2)*pi**2*t1)
				u_ana_2[i,j] += (-1)**(n+m-2)/(float(n*m))*sin(n*pi*x[i])*sin(m*pi*y[j])*exp(-(n**2 + m**2)*pi**2*t2)
		u_ana_1[i,j] = x[i]*y[j] - 4/pi**2*u_ana_1[i,j]
		u_ana_2[i,j] = x[i]*y[j] - 4/pi**2*u_ana_2[i,j]


x,y = meshgrid(x,y)


fig1 = figure()
ax1 = fig1.gca(projection = '3d')
surf = ax1.plot_surface(x,y,u1, rstride = 5, cstride = 5, cmap=cm.coolwarm,linewidth = 0,antialiased = False)
ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'$y$')
ax1.set_zlabel(r'$u(x,y)$')
fig1.savefig('Output/Diffusion_2d_t1_us.png')
#show()

fig2 = figure()
ax2 = fig2.gca(projection = '3d')
surf = ax2.plot_surface(x,y,u2, rstride = 5, cstride = 5, cmap=cm.coolwarm,linewidth = 0,antialiased = False)
ax2.set_xlabel(r'$x$')
ax2.set_ylabel(r'$y$')
ax2.set_zlabel(r'$u(x,y)$')
fig2.savefig('Output/Diffusion_2d_t2_us.png')
#show()


fig_d1 = figure()
ax_d1 = fig_d1.gca(projection = '3d')
surf = ax_d1.plot_surface(x,y,abs(u1-u_ana_1),rstride = 5, cstride = 5, cmap=cm.coolwarm,linewidth = 0,antialiased = False)
ax_d1.set_xlabel(r'$x$')
ax_d1.set_ylabel(r'$y$')
#ax_d1.set_zlabel(r'$\Delta u$')
fig_d1.savefig('Output/Delta_2d_t1_us.png')
#show()

fig_d2 = figure()
ax_d2 = fig_d2.gca(projection = '3d')
surf = ax_d2.plot_surface(x,y,abs(u2-u_ana_2),rstride = 5, cstride = 5, cmap=cm.coolwarm,linewidth = 0,antialiased = False)
ax_d2.set_xlabel(r'$x$')
ax_d2.set_ylabel(r'$y$')
#ax_d2.set_zlabel(r'$\Delta u$')
fig_d2.savefig('Output/Delta_2d_t2_us.png')
#show()

