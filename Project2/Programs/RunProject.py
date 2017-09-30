import os, sys

#omega_r = ['0.01', '0.5', '1', '5']
#rho_N = {'0.01':'0.02', '0.5':'4', '1':'10', '5':'60'}
value_combs = {'0.01':50, '0.5':10.0, '1':10.0, '5':2.0}

for i in value_combs.keys():
	#rho_N = str(5.0/float(i)) 
	os.system("./Project2.x 2pNoInt 200 "+str(value_combs[i])+" "+i)
	os.system("./Project2.x 2pCoulomb 200 "+str(value_combs[i])+" "+i)


