from math import * 
from matplotlib.pyplot import * 
from numpy import * 

omega_r = ['0.01', '0.5', '1', '5']

NoInt = {}
Coulomb = {}

for i in omega_r: 
	infile1 = open('Output/Eigenvectors_2pNoInt_omega'+i+'.txt', 'r')
	infile2 = open('Output/Eigenvectors_2pCoulomb_omega'+i+'.txt', 'r')
	
	NoInt[i] = []
	Coulomb[i] = []

	for line in infile1:  
		words = line.split();
		NoInt[i].append(float(words[0])) 

	for line in infile2:  
		words = line.split();
		Coulomb[i].append(float(words[0])) 	

	NoInt[i] = array(NoInt[i])
	Coulomb[i] = array(Coulomb[i])

	infile1.close() 
	infile2.close()



n = len(NoInt['1']) 

h = 10.0/n 
rho = [(i+1)*h for i in range(0,n)]
rho = array(rho)

x_range = {'0.01':12, '0.5':6, '1':5, '5':2}

for i in omega_r: 
	figure()
	plot(rho, NoInt[i]**2, 'r' )
	plot(rho, Coulomb[i]**2, 'b') 
	legend(['No interaction', 'Coulomb repulsion'])
	xlabel(r'$\rho$')
	ylabel(r'$u(\rho)^2$')
	xlim([0,x_range[i]])
	savefig('Output/TwoParticle_Eigenvectors_omega'+i+'.png')





