from math import * 
from matplotlib.pyplot import * 
from numpy import * 

font = {'size':14}
matplotlib.rc('font', **font) 

infile = open('Output/Eigenvectors_1pHO.txt', 'r')

R0 = []
R1 = []
R2 = []

n = 0; 

for line in infile:
	n += 1 
	words = line.split() 
	R0.append(float(words[0]))
	R1.append(float(words[1]))
	R2.append(float(words[2]))

R0 = array(R0)
R1 = array(R1)
R2 = array(R2)

infile.close()

h = 5.0/n 
rho = [(i+1)*h for i in range(0,n)]
rho = array(rho)


plot(rho, R0**2, 'b')
plot(rho, R1**2, 'r')
plot(rho, R2**2, 'g')
legend([r'$\lambda =3$',r'$\lambda =7$',r'$\lambda =11$'])
xlabel(r'$\rho$')
ylabel(r'$|u(\rho)|^2$')
xlim([0,5])
savefig('Output/Eigenvectors_1pHO.png')
#show()


