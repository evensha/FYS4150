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
		double tau, t; 
		tau = (A(l,l) - A(k,k))/(2*A(k,l)); 
		
		if(tau > 0){ t = - tau + sqrt(1.0 + tau*tau ); } // ????
		else{ t =  - tau - sqrt(1.0 + tau*tau ); } 
		c = 1.0/sqrt(1.0 + t*t);
		s = c*t; 
	}
	else{ c = 1.0; s = 0.0; } 

	double a_kk, a_ll, a_ik, a_il, r_ik, r_il; 
	a_kk = A(k,k); 
	a_ll = A(l,l);
	A(k,k) = a_kk*c*c - 2.0*A(k,l)*c*s + a_ll*s*s; 
	A(l,l) = a_kk*s*s + 2.0*A(k,l)*c*s + a_ll*c*c; 
	A(k,l) = 0.0; 
	A(l,k) = 0.0; 

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

	// Do the iterations 

	while( max_offdiag  > tolerance && iterations <= max_iterations ){  // iteration loop  
		Jacobi_rotation(A, R, p, q, n);
		max_offdiag = offdiag(A, &p, &q, n);
		iterations++;
	}

	cout << "Number of iterations: " << iterations << endl; 

	// Put eigenvalues in a vector

	for(int i = 0; i < n; i++ ){
			lambda(i) = A(i,i) ; 
	}
	
	return;  
}


void Jacobi_tests(mat A, int n){   // To be called if you want to test the algorithm before running it!    

	cout << "-------------------------------------" << endl; 
	cout << "OUTPUT FROM TEST FUNCTION:" << endl; 

	// Max off-diagonal element  

 	int p, q; 
	double test_max = offdiag(A, &p, &q, n); 

	cout << "Max off-diagonal element: " << A(p,q) << "   k,l=" << p << "," << q  << endl; 

	// Eigenvalues 

	mat R = zeros<mat>(n,n);
	vec l = zeros<vec>(n); 
	do_Jacobi(A,R,l,n);  

	l = sort(l); 

	cout << "Eigenvalues: "; 
	for(int i = 0; i<n; i++){
		cout << l(i) << " "; 
	}
	cout << endl; 	

	// Dot product of eigenvectors 

	vec R1(n); vec R2(n); 
	double c = 0;
	for(int i = 0; i<n; i++){  // pick out two first eigenvectors 
		R1(i) = R(i,0);  
		R2(i) = R(i,1); 
		c += R1(i)*R2(i); // dot product
		//cout << R1(i) << " , " << R2(i) << endl; 
	} 
 
	cout << "Dot product: " << c << endl; 
	cout << "-------------------------------------" << endl; 

	return; 
} 





