#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <iomanip>
//#include <armadillo>

using namespace std;

ofstream ofile; 

void Jacobi_rotation( double **A, double **D, int k, int l, int n ); // function that does the rotation 
void offdiag( double **A, int *p, int *q, int n ); // function that check value of the off-diagonal elements 

int main(int argc, char *argv[]){ 
	double tolerance = 1E-10; 
	int iterations = 0; 
	int n = atoi(argv[1]); 
	int max_iterations = n*n*n; 

	// Define vectors and matrces 
	double 

	// Defining eigenvector matrix 
	double D[n+1][n+1]; 
	for(int i = 0; i < n+1; i++ ){
		for(int j = 0; j < n+1; j++ ){
			if(i == j){D[i][j] = 1.0;} 
			else{D[i][j] = 0.0;}  
		}
	} 


/*
	while( max_offdiag > tolerance && iterations <= max_iterations ){  // iteration loop -> reduce off-diagonal elements to below the tolerance  
		int p, q; 
		offdiag(A, &p, &q, n);
		Jacobi_rotate(A, D, p, q, n);
		iterations++;
	}
*/

	return 0; 
} 


void offdiag( double **A, int *p, int *q, int n ){ 
	

} 
/*
void Jacobi_rotation( double **A, double **D, int k, int l, int n ){ 

}
*/


