import os, sys

omega_r = ['0.01', '0.5', '1', '5']

for i in omega_r: 
	os.system("./Project2.x 2pNoInt 400 10 "+i)
	os.system("./Project2.x 2pCoulomb 400 10 "+i)


