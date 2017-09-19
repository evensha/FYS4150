#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace std;
using namespace arma; 

ofstream ofile; 

void Jacobi_rotation( mat A, mat R, int k, int l, int n ); // function that does the rotation 
void offdiag( mat A, int *p, int *q, int n ); // function that check value of the off-diagonal elements 

int main(int argc, char *argv[]){ 
	double tolerance = 1E-10; 
	int iterations = 0; 
	int n = atoi(argv[1]); 
	int max_iterations = n*n*n; 

	double rho_N = atof(argv[2]); // ??????
	double h = rho_N/( (double) n);

	// Define rho and A
 
	vec rho(n); 
	mat A(n, n); 

	for(int i = 0; i < n; i++ ){
		rho(i) = i*i*h*h;		
	}

	for(int i = 0; i < n; i++ ){
		for(int j = 0; j < n; j++ ){	
			if(i == j){A(i,j) = 2.0/(h*h) + rho(i);}
			else if(i == j-1 or i == j+1){A(i,j) = -1.0/(h*h);} 
			else{A(i,j) = 0.0;}

			cout << A(i,j) << "  " ; if(j == n-1){cout << endl;}
		}
	}

	//vec eigval = eig_sym(A); 
	//for(int i = 0; i<3; i++){cout << eigval(i) << endl;}

	// Defining eigenvector matrix 
	mat R(n, n); 
	for(int i = 0; i < n; i++ ){
		for(int j = 0; j < n; j++ ){
			if(i == j){R(i,j) = 1.0;} 
			else{R(i,j) = 0.0;}  
		}
	} 

	
	//Jacobi_rotation(A, R, 1, 1, n); 

/*
	while( max_offdiag > tolerance && iterations <= max_iterations ){  // iteration loop -> reduce off-diagonal elements to below the tolerance  
		int p, q; 
		offdiag(A, &p, &q, n);
		Jacobi_rotate(A, R, p, q, n);
		iterations++;
	}
*/

	return 0; 
} 


void offdiag( double **A, int *p, int *q, int n ){   
	// start here tomorrow !!!!!!!
} 


void Jacobi_rotation( mat A, mat R, int k, int l, int n ){  
	double s,c;

	if( A(k,l) != 0){ 	
		cout << "A(k,l) was different from 0!" << endl; 
		double tau, t; 
		tau = (A(l,l) - A(k,k))/(2*A(k,l)); 
		
		if(tau >= 0){ t = - tau + sqrt(1+ tau*tau); } // ????
		else{ t = - tau - sqrt(1+ tau*tau); } 

		c = 1.0/sqrt(1 + t*t); 
		s = c*t; 
	}
	else{ c = 1.0; s = 0.0; } 

	double a_kk, a_ll, a_ik, a_il, r_ik, r_il; 
	a_kk = A(k,k); 
	a_ll = A(l,l);
	A(k,k) = a_kk*c*c - 2*A(k,l)*c*s + a_ll*s*s; 
	A(l,l) = a_kk*c*c + 2*A(k,l)*c*s + a_ll*s*s; 
	A(k,l) = 0.0; 
	A(l,k) = 0.0; 

	for( int i = 0; i < n; i++ ){ 
		if( i != k && i != l ){ 
			a_ik = A(i,k); 
			a_il = A(i,l);
			A(i,k) = a_ik*c + a_il*s;
			A(k,i) = A(i,k);  
			A(i,l) = a_il*c + a_ik*s; 
			A(l,i) = A(i,l); 
		}

		r_ik = R(i,k); 
		r_il = R(i,l); 

		R(i,k) = c*r_ik - s*r_il; 
		R(i,l) = c*r_il + s*r_ik; 
		cout << R(i,k) << " , " << R(i,l) << endl; 

	}
	return;  
}



