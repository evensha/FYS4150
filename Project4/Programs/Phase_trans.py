from matplotlib.pyplot import *
from numpy import *
import sys

font = {'size':18}
matplotlib.rc('font', **font)

files = {'40':open('Output/Output_L_40_dt0.01.txt','r'), '60':open('Output/Output_L_60_dt0.01.txt','r'), 
'80':open('Output/Output_L_80_dt0.01.txt','r'), '100':open('Output/Output_L_100_dt0.01.txt','r')}

T = []

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
		if L == '40':
			T.append(float(words[0]))  
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

L_values = []
for i in files.keys(): 
	L_values.append(int(i))

L_values.sort()


figure()
for L in L_values:
	plot(T,C_v[str(L)])
legend(['L=40', 'L=60', 'L=80', 'L=100'], loc = 2)
xlabel('T')
ylabel(r'$C_V$')
xlim(2,2.3)
savefig("Output/CV_T.png", bbox_inches = 'tight')

figure()
for L in  L_values:
	plot(T,Chi[str(L)])
legend(['L=40', 'L=60', 'L=80', 'L=100'], loc = 2)
xlabel('T')
ylabel(r'$\chi$')
xlim(2,2.3)
savefig("Output/Chi_T.png", bbox_inches = 'tight')

figure()
for L in L_values:
	plot(T,E[str(L)])
legend(['L=40', 'L=60', 'L=80', 'L=100'], loc = 2)
xlabel('T')
ylabel('E')
xlim(2,2.3)
savefig("Output/E_T.png", bbox_inches = 'tight')

figure()
for L in L_values:
	plot(T,AbsM[str(L)])
legend(['L=40', 'L=60', 'L=80', 'L=100'], loc = 2)
xlabel('T')
ylabel('|M|')
xlim(2,2.3)
savefig("Output/M_T.png", bbox_inches = 'tight')


