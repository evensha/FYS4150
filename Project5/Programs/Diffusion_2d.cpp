#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <string>
#include <iomanip>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma; 

ofstream ofile; 

int main( int argc, char *argv[] ){

	int n = atoi(argv[1]); 
	double h = 1.0/((double) n +1.0); 
	double delta_t = h*h/2.0; 
	double t_f = 10/delta_t;
	double alpha = 0.5*delta_t/(h*h);  
	double t1 = 0.05/delta_t; double t2 = 0.5/delta_t; 
	double epsilon = 0.00001; 
	int iterations; 
	int max_iterations = 10000; 
	double diff = 0; 
	double delta_u = 0; 

	mat u = zeros<mat>(n+1,n+1); 
	mat unew = zeros<mat>(n+1,n+1); 
	mat u_temp = zeros<mat>(n+1, n+1); 
	mat u_prev = zeros<mat>(n+1, n+1); 
	mat u_t1; 

	//u(n,n) = 1.0 
	for( int i = 0; i <= n; i++ ){ // set up boundary conditions 
		u(n,i) = i/((double) n); 
		u(i,n) = i/((double) n); 
	}
	unew = u; 

	cout << h << endl; 
	cout << u(n,n) << endl; 
	cout << u(0,0) << endl; 

	for( int t = 1; t <= t_f; t++ ){	// loop over t
		iterations = 0;  
		//u_prev = u; 
		u = unew; 
		//while( (iterations <= max_iterations) && (diff > epsilon) ){	// while loop over iterations 
			//u_temp = u; diff = 0;  
			for( int i = 1; i < n; i++ ){	// loop over x 
				for( int j = 1; j < n; j++ ){	// loop over y 
					// update something 
					//delta_u = u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1); 
					//u(i,j) = 1.0/(1.0 + 4.0*alpha)*(alpha*delta_u + u_prev(i,j)); 
					//diff += fabs(u_temp(i,j)-u(i,j)); 
					unew(i,j) = u(i,j) + alpha*( u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4.0*u(i,j)); 
				}
			}	
			//iterations++;
			//diff /= pow((n),2.0);  			
		//}
		if(t == t1) u_t1 = u; 
	}
		
	ofile.open("Output/Diffusion_2d.txt"); 
	for( int i = 0; i <= n; i++ ){ 
		for( int j = 0; j <= n; j++ ){ 
			ofile << u(i,j) << setw(20); 
		}
		ofile << endl; 
	}
	ofile.close(); 

	return 0; 

}


