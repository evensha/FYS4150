#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <string>
#include <iomanip>
#include <armadillo>
#include "time.h"

#include "Jacobi_algorithm.h"

using namespace std;
using namespace arma; 

ofstream ofile; 


int main(int argc, char *argv[]){ 

	// First test matrix (eigenvalues should be 1 and 6)

	mat T1(2,2); 
	T1(0,0) = 5.0; T1(0,1) = -2.0; T1(1,0) = -2.0; T1(1,1) = 2.0; 
	Jacobi_tests(T1,2); // Check that the algorithm works as it should  

	// Second test matrix 

	mat T2 = ones<mat>(5,5); T2(1,2) = 5.0; 
	Jacobi_tests(T2,5);

	// Which problem should be considered? (1pHO, 2pNoInt, 2pCoulomb)

	string Prob = argv[1]; 

	if(!(Prob == "1pHO" or Prob == "2pNoInt" or Prob == "2pCoulomb")){
		cout << "UNKNOWN PROBLEM!" << endl; 
		exit(1); 
	} 

	// Define some constants

	int n = atoi(argv[2]); 

	double rho_0 = 0.0; 
	double rho_N = atof(argv[3]); // 10.0; 
	double h = (rho_N - rho_0)/( (double) n );
	double omega_r; 

	if( Prob == "1pHO" ){ omega_r = 1.0; } 
	else{ omega_r = atof(argv[4]); }   

	double diag = 2.0/(h*h); 
	double non_diag = -1.0/(h*h); 


	// Define necessary matrices and vectors 


	vec V = zeros<vec>(n); // The potential 

	for(int i = 0; i < n; i++ ){
		double rho_i = rho_0 + (i+1.0)*h;	
		if( Prob == "1pHO" or Prob == "2pNoInt" ){	V(i) = omega_r*omega_r*rho_i*rho_i; }  
		if( Prob == "2pCoulomb" ){ V(i) = omega_r*omega_r*rho_i*rho_i + 1.0/rho_i; }
	}

	mat  A = zeros<mat>(n, n);  // The A matrix (which we want to diagonalize)  

	A(0,0) = diag + V(0); 
	A(0,1) = non_diag; 

	for(int i = 1; i < n-1; i++ ){
		A(i,i) = diag + V(i); 
		A(i,i-1) = non_diag; 
		A(i,i+1) = non_diag; 					
	}

	A(n-1, n-2) = non_diag; 
	A(n-1, n-1) = diag + V(n-1); 

	mat R = zeros<mat>(n,n); 	 // Matrix that will contain the eigenvectors 
	vec lambda = zeros<vec>(n);  // Vector that will contain the eigenvalues 


	// Eigenvalues from Armadillo 
/*
	vec eigval = eig_sym(A);      
	cout << "Eigenvalues from Armadillo:" << endl; 
	for(int i = 0; i < n; i++){
		cout << eigval(i) << endl;
	}
*/


	// Run the Jacobi algorithm

	clock_t start, finish; 
	start = clock(); 

	do_Jacobi(A, R, lambda, n); 

	finish = clock(); 
	double time = (finish -start)/((double) CLOCKS_PER_SEC); 

	if(time > 60){ cout << "Time spent: " << time/60.0 << " min" << endl;} 
	else cout << "Time spent: " << time << " s" << endl;
	

	// Process the output from the Jacobi algorithm 

	map<double, vec> Eigen;   // make a map to pair up eigenvectors and eigenvalues

	for(int i = 0; i < n; i++){   // run over the matrix R to pick out eigenvectors 
		vec Eigenvector = zeros<vec>(n); 
		for(int j = 0; j < n; j++){
			Eigenvector(j) = R(j,i); 
		}  
		Eigen[lambda(i)] = Eigenvector;  // pair up eigenvalues and eigenvectors 
	}  

	lambda = sort(lambda);   // sort eigenvalues 

	cout << "Three lowest eigenvalues:" << endl; 
	cout << lambda(0) << endl; 
	cout << lambda(1) << endl; 
	cout << lambda(2) << endl; 
	cout << "-------------------------------------" << endl; 


	// Write eigenvectors to file 

	ostringstream os; 
	if( Prob == "1pHO" ){
		os << "Output/Eigenvectors_" << Prob << ".txt"; 
	}
	else{
		os << "Output/Eigenvectors_" << Prob << "_omega" << omega_r << ".txt"; 
	}

	string outfile = os.str(); 
	ofile.open(outfile.c_str()); 

	for(int i = 0; i < n; i++){
		ofile << Eigen[lambda(0)](i) << setw(20) << Eigen[lambda(1)](i) << setw(20) << Eigen[lambda(2)](i) << endl; 						
	}	 

	ofile.close(); 

	return 0; 
} 







