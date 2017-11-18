from matplotlib.pyplot import *
from numpy import *
import sys

#infile40 = open('Output/Output_L_40.txt','r')
#infile60 = open('Output/Output_L_60.txt','r')
#infile80 = open('Output/Output_L_80.txt','r')
#infile100 = open('Output/Output_L_100.txt','r')

files = {'40':open('Output/Output_L_40.txt','r'), '60':open('Output/Output_L_60.txt','r'), '80':open('Output/Output_L_80.txt','r'),
'100':open('Output/Output_L_100.txt','r')}
#files = {'80':open('Output/Output_L_80.txt','r')}
T = [2.0, 2.025, 2.05, 2.075, 2.1, 2.125, 2.15, 2.175, 2.2, 2.225, 2.25, 2.275, 2.3]

E = {}
C_v = {}
M = {}
Chi = {}
AbsM = {}

for L in files.keys(): 
	E[L] = []
	C_v[L] = []
	M[L] = []
	Chi[L] = []
	AbsM[L] = []
	for line in files[L]: 
		words = line.split()
		E[L].append(float(words[1])) 
		C_v[L].append(float(words[2]))
		M[L].append(float(words[3]))
		Chi[L].append(float(words[4]))
		AbsM[L].append(float(words[5]))			
	E[L] = array(E[L])
	C_v[L] = array(C_v[L])
	M[L] = array(M[L])
	Chi[L] = array(Chi[L])
	AbsM[L] = array(AbsM[L])


figure()
for L in files.keys():
	plot(T,C_v[L])
show()









