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

int main(int argc, char *argv[]){ 

	string method = argv[1]; 

	int n = atoi(argv[2]); 
	double t_0 = 0.0; 
	double t_n = 5.0; 
	double h = (t_0 - t_n)/((double) n); 
	
	vec x = zeros<vec>(n); 
	vec y = zeros<vec>(n); 
	vec v_x = zeros<vec>(n); 
	vec v_y = zeros<vec>(n); 
	vec a_x = zeros<vec>(n); 
	vec a_y = zeros<vec>(n); 
	
	// Initial conditions (from NASA) 

	double x_0 = 9.734E-01; double y_0 = 2.4198E-01; 
	double vx_0 = -4.349E-03*365; double vy_0 = 1.665E-02*365;  

	x(0) = x_0; y(0) = y_0; v_x(0) = vx_0; v_y(0) = vy_0; 

	// Forward Euler 

	if(method == "FE"){
		for(int i = 0; i<n-1; i++){
			double r = sqrt(x(i)*x(i) + y(i)*y(i)); 
			a_x(i) = - 4.0*M_PI*M_PI*x(i)/(r*r*r); 
			a_y(i) = - 4.0*M_PI*M_PI*y(i)/(r*r*r);

			x(i+1) = x(i) + h*v_x(i); 
			v_x(i+1) = v_x(i) + h*a_x(i); 
		
			y(i+1) = y(i) + h*v_y(i); 
			v_y(i+1) = v_y(i) + h*a_y(i);   
		}
	}

	// Velocity Verlet  

	if(method == "VV"){

		double r_0 = sqrt(x_0*x_0 + y_0*y_0); 
		a_x(0) = - 4.0*M_PI*M_PI*x_0/(r_0*r_0*r_0); 
		a_y(0) = - 4.0*M_PI*M_PI*y_0/(r_0*r_0*r_0); 

		for(int i = 0; i<n-1; i++){

			x(i+1) = x(i) + h*v_x(i) + h*h/2.0*a_x(i); 
			y(i+1) = y(i) + h*v_y(i) + h*h/2.0*a_y(i); 

			double r = sqrt(x(i+1)*x(i+1) + y(i+1)*y(i+1)); 
			
			a_x(i+1) = - 4.0*M_PI*M_PI*x(i+1)/(r*r*r); 	
			a_y(i+1) = - 4.0*M_PI*M_PI*y(i+1)/(r*r*r);		

			v_x(i+1) = v_x(i) + h/2.0*(a_x(i+1) + a_x(i));
			v_y(i+1) = v_y(i) + h/2.0*(a_y(i+1) + a_y(i));

		}
	}
		
	// Write to file 

	ostringstream os; 
	os << "Output/EarthSun_"<< method  <<".txt"; 

	string outfile = os.str(); 
	ofile.open(outfile.c_str()); 

	for(int i = 0; i < n; i++){
		ofile << x(i) << setw(20) << y(i) << setw(20) << v_x(i) << setw(20) << v_y(i) << endl; 						
	}	 

	ofile.close();

} 
