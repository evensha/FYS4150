#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma; 

ofstream ofile; 


double offdiag( mat A, int *p, int *q, int n ); // function that check value of the off-diagonal elements 
void Jacobi_rotation( mat& A, mat& R, int k, int l, int n ); // function that does the rotation 


int main(int argc, char *argv[]){ 
	double tolerance = 1.0E-10; 
	int iterations = 0; 
	int n = atoi(argv[1]); 
	int max_iterations = n*n*n; 

	double rho_0 = 0.0; 
	double rho_N = 10.0; // sqrt(n*n + 1.0 ); // sqrt(10.0); // atof(argv[2]); // ??????
	double h = (rho_N - rho_0)/( (double) n );

	double diag = 2.0/(h*h); 
	double non_diag = -1.0/(h*h); 

	// Define rho, V and A

	vec rho(n); vec V(n); 
	mat  A = zeros<mat>(n, n); 

	for(int i = 0; i < n; i++ ){
		rho(i) = rho_0 + (i+1.0)*h;		
		V(i) = rho(i)*rho(i); 
		//cout <<  rho(i-1) << ", " << endl;
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

	//cout << "------------" << endl; 
/*
	for(int i = 0; i < n; i++ ){
		for(int j = 0; j < n; j++ ){
			cout <<	A(i,j) << "  "; if(j==n-1){cout << endl;}	
		}
	}
*/
	//mat A(2,2); 
	//A(0,0) = 5.0; A(0,1) = -2.0; A(1,0) = -2.0; A(1,1) = 2.0; 


	// Eigenvalues from Armadillo 
/*
	vec eigval = eig_sym(A); 
	cout << "Eigenvalues from Armadillo:" << endl; 
	for(int i = 0; i < n; i++){
		cout << eigval(i) << endl;
	}
	cout << "--------------" << endl; 
*/
	
	// Define eigenvector matrix
 
	mat R = zeros<mat>(n, n); 
	for(int i = 0; i < n; i++ ){
		R(i,i) = 1.0; 
	} 


	// Run the algorithm  

	int p, q; 
	double max_offdiag = offdiag(A, &p, &q, n); 

	clock_t start, finish; 
	start = clock(); 

	while( max_offdiag  > tolerance && iterations <= max_iterations ){  // iteration loop  
		//cout << "p,q = " << p << " , " << q << endl;
		Jacobi_rotation(A, R, p, q, n);
		max_offdiag = offdiag(A, &p, &q, n);
		iterations++;
		//cout << max_offdiag << endl; 
	}

	finish = clock(); 
	double time = (finish -start)/((double) CLOCKS_PER_SEC); 
	cout << "Time spent: " << time << endl; 
	cout << "Number of iterations: " << iterations << endl;
 

	// Eigenvalues

	vec lambda(n); 
	cout << "Three lowest eigenvalues:" << endl; 
	for(int i = 0; i < n; i++ ){
			lambda(i) = A(i,i) ; 
	}

	lambda = sort(lambda); 

	cout << lambda(0) << endl; 
	cout << lambda(1) << endl; 
	cout << lambda(2) << endl; 

	return 0; 
} 


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



