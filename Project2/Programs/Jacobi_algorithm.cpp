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

ofstream ofile1, ofile2; 

//--------------------------------------
// Maximal off-diagonal matrix element  
//--------------------------------------

double offdiag( mat A, int *p, int *q, int n ){   
	double max = 0.0; 
	for(int i = 0; i<n; i++){
		for(int j = i+1; j<n; j++){
			double a_ij = fabs(A(i,j)); 
			if(a_ij > max){ 
				max = a_ij; *p = i; *q = j; 	
			}
		} 
	}
	return max; 
} 

//------------------
// Jacobi rotation  
//------------------

void Jacobi_rotation( mat& A, mat& R, int k, int l, int n ){  
	double s,c;

	if( A(k,l) != 0){ 	
		//cout  << "A(k,l) = " << A(k,l) << endl; 
		double tau, t; 
		tau = (A(l,l) - A(k,k))/(2*A(k,l)); 
		
		if(tau > 0){ t =  1.0/( tau + sqrt(1.0 + tau*tau )); } // ????
		else{ t = -1.0/( - tau + sqrt(1.0 + tau*tau )); } 
		//cout << "tau,t = " << tau << " , " << t << endl; 
		c = 1.0/sqrt(1.0 + t*t);
		s = c*t; 
		//cout << "c,s = " << c << "  " << s << endl; 
	}
	else{ c = 1.0; s = 0.0; } 

	double a_kk, a_ll, a_ik, a_il, r_ik, r_il; 
	a_kk = A(k,k); 
	a_ll = A(l,l);
	/*
	cout << "-----------------" << endl;	
	cout << "Before rotation:" << endl;
	cout << "l,k = " << l << "," << k << endl;
	cout << "s,c = " << s << "," << c << endl;
	cout << "A(k,k) = " << A(k,k) << endl; 
	cout << "A(l,l) = " << A(l,l) << endl; 
	cout << "-----------------" << endl;
	*/
	A(k,k) = a_kk*c*c - 2.0*A(k,l)*c*s + a_ll*s*s; 
	A(l,l) = a_kk*s*s + 2.0*A(k,l)*c*s + a_ll*c*c; 
	A(k,l) = 0.0; 
	A(l,k) = 0.0; 
	/*
	cout << "-----------------" << endl;
	cout << "After rotation:" << endl; 
	cout << "A(k,k) = " << A(k,k) << endl; 
	cout << "A(l,l) = " << A(l,l) << endl; 
	cout << "-----------------" << endl; 
	*/
	for( int i = 0; i < n; i++ ){ 
		if( i != k && i != l ){ 
			a_ik = A(i,k); 
			a_il = A(i,l);
			A(i,k) = a_ik*c - a_il*s;
			A(k,i) = A(i,k);  
			A(i,l) = a_il*c + a_ik*s; 
			A(l,i) = A(i,l); 
		}

		r_ik = R(i,k); 
		r_il = R(i,l); 

		R(i,k) = c*r_ik - s*r_il; 
		R(i,l) = c*r_il + s*r_ik; 
		//cout << R(i,k) << " , " << R(i,l) << endl; 

	}
	return;  
}

//------------------------
// Running the algorithm 
//------------------------

void do_Jacobi(mat& A, mat& R, vec& lambda, int n){

	// Define the eigenvector matrix
 
	//mat R = zeros<mat>(n, n); 
	for(int i = 0; i < n; i++ ){
		R(i,i) = 1.0; 
	} 

	// Define some necessary quantities 

	double tolerance = 1.0E-10; 
	int max_iterations = n*n*n; 
	int iterations = 0; 

	int p, q; 
	double max_offdiag = offdiag(A, &p, &q, n);  // find the largest element of A (i.e. starting point of the iteration loop) 

	clock_t start, finish; 
	start = clock(); 

	// Do the iterations 

	while( max_offdiag  > tolerance && iterations <= max_iterations ){  // iteration loop  
		Jacobi_rotation(A, R, p, q, n);
		max_offdiag = offdiag(A, &p, &q, n);
		iterations++;
	}

	// Put eigenvalues in a (sorted) vector

	for(int i = 0; i < n; i++ ){
			lambda(i) = A(i,i) ; 
	}

	lambda = sort(lambda); 

	// Write results to file 
/*
	ofile1.open("Eigenvalues.txt");
	ofile2.open("Eigenvectors.txt"); 

	for(int i = 0; i < n; i++){ 
		ofile1 << "Lambda_" << i << " = " << lambda(i) << endl;  
		for(int j = 0; j < n; j++){ 
			ofile2 << R(i,j) << "  " ; if(j == n-1){ofile2 << endl;}
		}
	}

	ofile1.close(); 
	ofile2.close(); 
*/
	// Finish... 

	finish = clock(); 
	double time = (finish -start)/((double) CLOCKS_PER_SEC); 
	if( n > 2 ){ 
		cout << "Time spent: " << time << endl; 
		cout << "Number of iterations: " << iterations << endl;
	}

	return;  
}


void Jacobi_tests(){   // To be called if you want to test the algorithm before running it!  

	// Max off-diagonal element of a 2x2 matrix 

	mat T2(2,2); 
	T2(0,0) = 1.0; T2(0,1) = 10.0; T2(1,0) = 5.0; T2(1,1) = 1.0;
 	int p, q; 
	double test_max = offdiag(T2, &p, &q, 2); 

	if(test_max != 10.0 or p != 0 or q != 1){  
		cout << "DID NOT PASS MAX OFF-DIAGONAL ELEMENT TEST!" << endl; exit(1);
	}	

	// Eigenvalues of 2x2 matrix 

	mat T1(2,2); 
	T1(0,0) = 5.0; T1(0,1) = -2.0; T1(1,0) = -2.0; T1(1,1) = 2.0; // setting up a matrix with known eigenvalues 1 and 6 
	mat S = zeros<mat>(2,2);
	vec l = zeros<vec>(2); 
	do_Jacobi(T1,S,l,2);  

	if(l(0) != 1.0 or l(1) != 6.0){   
		cout << "DID NOT PASS EIGENVALUETEST!" << endl; exit(1); 
	}


	cout << "---------------------" << endl; 
	cout << "ALL TESTS PASSED! :-D" << endl; 
	cout << "---------------------" << endl; 	


	return; 
} 





