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

void Jacobi_tests(); 

int main(int argc, char *argv[]){ 

	// Call test function to check thing works as they should 

	Jacobi_tests(); 

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

	if( Prob == "1pHO" ){ omega_r = 1.0; }  // omega_r set to 1 when considering one-particle case
	else{ omega_r = atof(argv[4]); }   

	double diag = 2.0/(h*h); // pre-define the constant part of the diagonal elements 
	double non_diag = -1.0/(h*h); // pre-define the off-diagonal elements 


	// Define necessary matrices and vectors 

	vec V = zeros<vec>(n); // The potential 

	for(int i = 0; i < n; i++ ){
		double rho_i = rho_0 + (i+1.0)*h;	 // rho_0 is known so the first element of rho should be rho_1, hence the (i+1)*h
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


	// Solving with Armadillo 

	clock_t arma_start, arma_finish; 
	arma_start = clock(); 
	vec eigval = eig_sym(A);      
	arma_finish = clock(); 
	double arma_time = (arma_finish - arma_start)/((double) CLOCKS_PER_SEC); 
	cout << "Time with Armadillo: " << arma_time << " s" << endl; 
	/*
	cout << "Eigenvalues from Armadillo:" << endl; 
	for(int i = 0; i < 3; i++){
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


void Jacobi_tests(){   // To be called if you want to test the algorithm before running it!    

	// Test eigenvalues 

	mat T1(2,2); mat R1 = zeros<mat>(2,2); vec lambda1 = zeros<vec>(2); 
	T1(0,0) = 5.0; T1(0,1) = -2.0; T1(1,0) = -2.0; T1(1,1) = 2.0;    // matrix with eigenvalues 1 and 6

	do_Jacobi(T1, R1, lambda1, 2); 

	lambda1 = sort(lambda1); 	
	if(!(lambda1(0) == 1 && lambda1(1) == 6)){ cout << "EIGENVALUE TEST NOT PASSED!" << endl; exit(1);}  


	// Test max off-diagonal element 

	mat T2 = ones<mat>(5,5); T2(1,2) = 5.0; T2(4,3) = 3.0; 
 	int p, q; 
	double test_max = offdiag(T2, &p, &q, 5); 

	if(test_max != 5.0){cout << "MAX OFF-DIAG ELEMENT TEST NOT PASSED!" << endl; exit(1);} 


	// Test dot product 

	mat R2 = zeros<mat>(5,5); vec lambda2 = zeros<vec>(5);  
	do_Jacobi(T2,R2,lambda2,5); // using the same matrix as above    	

	vec V1(5); vec V2(5); 
	double c1 = 0;
	double c2 = 0; 

	for(int i = 0; i<5; i++){  // pick out two first eigenvectors 
		V1(i) = R2(i,0);  
		V2(i) = R2(i,1); 
		c1 += V1(i)*V2(i); // dot product
		c2 += V1(i)*V1(i); 
	} 

	cout << c1 << endl; 
	cout << c2 << endl; 

	if( c1 > 1E-10 && c2 != 1.0 ){cout << "DOT PRODUCT TEST NOT PASSED!" << endl; exit(1);}  

	cout << "ALL UNIT TESTS PASSED! :-)" << endl; 
	cout << "--------------------------" << endl; 
 
	return; 
} 





