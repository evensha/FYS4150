import os, sys

value_combs = {'0.01':50, '0.5':6.0, '1':5.0, '5':2.0}
N = {'0.01':500, '0.5':400, '1':400, '5':200}

for i in value_combs.keys():
	print 'Running code with omega_r = %s' %i
	#rho_N = str(5.0/float(i)) 
	os.system("./Project2.x 2pNoInt "+str(N[i])+" "+str(value_combs[i])+" "+i)
	os.system("./Project2.x 2pCoulomb "+str(N[i])+" "+str(value_combs[i])+" "+i)


