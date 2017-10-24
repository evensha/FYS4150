#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <iomanip>
#include <string>
//#include <armadillo> 

using namespace std;
//using namespace arma; 

ofstream ofile; 

inline double f(double x){return 100.0*exp(-10.0*x);}
inline double exact(double x){return 1.0 - (1-exp(-10))*x -exp(-10*x);} 

int main(int argc, char* argv[]){

	string filename; 
	int n; 
	filename = argv[1];
	n = atoi(argv[2]); 

	double h = 1.0/(n+1.0); 
	double h2 = h*h;

	double *a = new double [n+2]; 
	double *b = new double [n+2]; 
	double *c = new double [n+2]; 

	double *x = new double [n+2]; 
	double *v = new double [n+2]; //solution 
	
	for(int i = 0; i < n+2; i++){ x[i] = i*h; }
	    				 		
	
	ofile.open(filename); 

	ofile.close(); 

	delete a; delete b; delete c; delete x; delete v; 
	return 0; 
}  
