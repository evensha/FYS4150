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

void Forward_Euler(double dx, double dt, int t_f, int t1, int t2, int n); 
void Backward_Euler(double dx, double dt, int t_f, int t1, int t2, int n); 
void Crank_Nicolson(double dx, double dt, int t_f, int t1, int t2, int n); 


int main( int argc, char *argv[] ){

	int n =  atoi(argv[1]); 
	double h = 1.0/((double) n +1.0); 
	double delta_t = h*h/2.0; 
	double t_f = 10/delta_t; 
	double t1 = 0.05/delta_t; double t2 = 0.5/delta_t; 
	cout << "------------" << endl; 
	cout << "n = " << n << endl; 
	cout << "h = " << h << endl;  
	cout << "dt = " << delta_t << endl;
	cout << "t_f = " << t_f << endl;
	cout << "t_1 = " << t1 << endl; 
	cout << "t_2 = " << t2 << endl;
	cout << "------------" << endl;        
	Forward_Euler(h, delta_t, t_f, t1, t2, n); 	
	Backward_Euler(h, delta_t, t_f, t1, t2, n); 
	Crank_Nicolson(h, delta_t, t_f, t1, t2, n); 

	return 0; 

}


void Forward_Euler(double dx, double dt, int t_f, int t1, int t2, int n){

	vec u = zeros<vec>(n+1);  				
	vec unew = zeros<vec>(n+1); 

	vec u_t1 = zeros<vec>(n+1); 
	vec u_t2 = zeros<vec>(n+1); 

	double alpha = dt/(dx*dx); 

	ofile.open("Output/Output_FE.txt"); 
	
	// Boundary conditions 
	u(0) = unew(0) = 0.0; 
	u(n) = unew(n) = 1.0; 

	// Initial condition 
	for( int i = 1; i<n; i++ ){ // this is not really necessary since vectors are zero anyway 
		u(i) = 0.0; 
		unew(i) = 0.0; 
	}
	//for( int j = 0; j<= n; j++ ) cout << j*dx << endl; 
	// Forward-Euler solver 
	for( int t = 1; t <= t_f; t++ ){
		for( int i = 1; i<n; i++ ){
			
			unew(i) = alpha*u(i-1) + (1.0 - 2.0*alpha)*u(i) + alpha*u(i+1); 	
			//u(i) = unew(i); 

		}
		
		if( t == t1 ) u_t1 = u; // keep results at t1 and t2
		if( t == t2 ) u_t2 = u; 
		u = unew; 
	}

	for( int i = 0; i<=n; i++ ){ // write u to file for t1 and t2 
		ofile << u_t1(i) << setw(20) << u_t2(i) << endl; 
	}

	ofile.close(); 
	return; 

}


void Backward_Euler(double dx, double dt, int t_f, int t1, int t2, int n){

	vec u, unew, d, b; 
	u = unew = zeros<vec>(n+1);
	d = zeros<vec>(n+1); 
	b = zeros<vec>(n+1);  	
	double alpha = dt/(dx*dx); 	
	vec u_t1 = zeros<vec>(n+1); 
	vec u_t2 = zeros<vec>(n+1); 

	ofile.open("Output/Output_BE.txt");

	// Boundary conditions 
	u(0) = unew(0) = 0.0; 
	u(n) = unew(n) = 1.0; 

	// Initial condition 
	for( int i = 1; i<n; i++ ){ // this is not really necessary since vectors are zero anyway 
		u(i) = 0.0; 
		unew(i) = 0.0; 
	}

	for( int t = 1; t<= t_f; t++ ){	

		for( int i = 1; i<n; i++ )  d(i) = 1.0+2.0*alpha; 
		for( int i = 1; i<n; i++ )	b(i) = -alpha;  

		// Forward substitution 
		for( int i = 2; i<n; i++ ){

			d(i) = d(i) - b(i-1)*b(i-1)/d(i-1); 
			u(i) = u(i) - b(i-1)*u(i-1)/d(i-1);

		}

		// Backward substitution 
		for( int i = n; i > 1; i-- ){
			unew(i-1) = ( u(i-1) - b(i-1)*unew(i) )/d(i-1); 
		}

		if( t == t1 ) u_t1 = unew; // keep results at t1 and t2
		if( t == t2 ) u_t2 = unew; 
		u = unew; 
	}

	for( int i = 0; i<=n; i++ ){ // write u to file for t1 and t2 
		ofile << u_t1(i) << setw(20) << u_t2(i) << endl; 
	}

	ofile.close(); 
	return; 

}


void Crank_Nicolson(double dx, double dt, int t_f, int t1, int t2, int n){

	vec u, unew, r, d, b; 
	u = unew = r = zeros<vec>(n+1);
	d = zeros<vec>(n+1); 
	b = zeros<vec>(n+1);  	
	double alpha = dt/(dx*dx); // (dt+dt/2.0)/(dx*dx); 	
	vec u_t1 =  zeros<vec>(n+1); 
	vec u_t2 = zeros<vec>(n+1); 

	ofile.open("Output/Output_CN.txt");

	// Boundary conditions 
	u(0) = unew(0) = r(0) = 0.0; 
	u(n) = unew(n) = r(n) = 1.0; 

	// Initial condition 
	for( int i = 1; i<n; i++ ){ // this is not really necessary since vectors are zero anyway 
		u(i) = 0.0; 
		unew(i) = 0.0; 
		r(i) = 0.0; 
	}
			
	for( int t = 1; t<= t_f; t++ ){

		for( int i = 1; i<n; i++ ){
			r(i) = alpha*u(i-1) + (2.0-2.0*alpha)*u(i) + alpha*u(i+1); 
		}

		for( int i = 1; i<n; i++ )  d(i) = 2.0+2.0*alpha; 
		for( int i = 1; i<n; i++ )	b(i) = -alpha;  

		// Forward substitution 
		for( int i = 2; i<n; i++ ){
			d(i) = d(i) - b(i-1)*b(i-1)/d(i-1); 
			r(i) = r(i) - b(i-1)*r(i-1)/d(i-1);
		}

		// Backward substitution 
		for( int i = n; i > 1; i-- ){
			unew(i-1) = ( r(i-1) - b(i-1)*unew(i) )/d(i-1); 
		}

		if( t == t1 ) u_t1 = unew; // keep results at t1 and t2
		if( t == t2 ) u_t2 = unew; 
		u = unew; 
	}

	for( int i = 0; i<=n; i++ ){ // write u to file for t1 and t2 
		ofile << u_t1(i) << setw(20) << u_t2(i) << endl; 
	}

	ofile.close(); 
	return; 

}



