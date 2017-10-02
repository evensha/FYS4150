from math import * 
from matplotlib.pyplot import * 
from numpy import * 

#omega_r = ['0.01', '0.5', '1', '5']
#omega_r = ['5']
value_combs = {'0.01':50.0, '0.5':10.0, '1':10.0, '5':2.0}


NoInt = {}
Coulomb = {}

for i in value_combs.keys(): 
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


	n = len(NoInt[i]) 


	rho_N = value_combs[i] 
	h = rho_N/n; 
	rho = [(j+1)*h for j in range(0,n)]
	rho = array(rho)


	figure()
	plot(rho, NoInt[i]**2, 'r' )
	plot(rho, Coulomb[i]**2, 'b') 
	legend(['No interaction', 'Coulomb repulsion'])
	xlabel(r'$\rho$')
	ylabel(r'$u(\rho)^2$')
	xlim([0, rho_N ])
	savefig('Output/TwoParticle_Eigenvectors_omega'+i+'.png')





