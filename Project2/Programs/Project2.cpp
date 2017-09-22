#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "time.h"

#include "Jacobi_algorithm.h"

using namespace std;
using namespace arma; 

//ofstream ofile; 


int main(int argc, char *argv[]){ 

	Jacobi_tests(); // Check that the algorithm works as it should  

	// Define some constants

	int n = atoi(argv[1]); 

	double rho_0 = 0.0; 
	double rho_N = 10.0; 
	double h = (rho_N - rho_0)/( (double) n );

	double diag = 2.0/(h*h); 
	double non_diag = -1.0/(h*h); 

	// Define rho, V and A

	vec rho(n); vec V(n); 
	mat  A = zeros<mat>(n, n); 

	for(int i = 0; i < n; i++ ){
		rho(i) = rho_0 + (i+1.0)*h;		
		V(i) = rho(i)*rho(i); 
	}

	A(0,0) = diag + V(0); 
	A(0,1) = non_diag; 

	for(int i = 1; i < n-1; i++ ){
		A(i,i) = diag + V(i); 
		A(i,i-1) = non_diag; 
		A(i,i+1) = non_diag; 					
	}

	A(n-1, n-2) = non_diag; 
	A(n-1, n-1) = diag + V(n-1); 


	// Eigenvalues from Armadillo 
/*
	vec eigval = eig_sym(A); 
	cout << "Eigenvalues from Armadillo:" << endl; 
	for(int i = 0; i < n; i++){
		cout << eigval(i) << endl;
	}
	cout << "--------------" << endl; 
*/
	

	
	mat R = zeros<mat>(n,n); 	
	vec lambda = zeros<vec>(n); 
	do_Jacobi(A, R, lambda, n); 


	//vec lambda(n); 

	cout << "Three lowest eigenvalues:" << endl; 

	cout << lambda(0) << endl; 
	cout << lambda(1) << endl; 
	cout << lambda(2) << endl; 

	return 0; 
} 







