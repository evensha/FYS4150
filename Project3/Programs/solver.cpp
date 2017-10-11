#include "solver.h"
#include "planet.h"
#include <iostream> 
#include <cmath>
#include "time.h"
#include <cstdlib> 
#include <fstream>
#include <iomanip>


using namespace std; 

ofstream ofile; 

solver::solver(double r){

	total_planets = 0; 
	mass = 0; 
	radius = r; 
	G = 4*M_PI*M_PI; 

}
	

void solver::addPlanet(planet newplanet){ 

	total_planets += 1; 
	all_planets.push_back(newplanet); 
	mass += newplanet.mass; 

} 


void solver::ForwardEuler(int integration_points, double time){

	// Constants & vectors 

	int n = integration_points; 	
	double h = time/((double) n);

	//double x[n]; double y[n]; double v_x[n]; double v_y[n]; double a_x[n]; double a_y[n]; 

	for(int i = 0; i<total_planets; i++){	 
		ofile.open("Output/EarthSun_FE.txt"); 

		planet &Planet = all_planets[i]; 

		double a_x; double a_y; 

		// Solve with Forward-Euler algo 

		for(int i = 0; i<n-1; i++){
			double r = sqrt(Planet.position[0]*Planet.position[0] + Planet.position[1]*Planet.position[1]); 
			a_x = - 4.0*M_PI*M_PI*Planet.position[0]/(r*r*r); 
			a_y = - 4.0*M_PI*M_PI*Planet.position[1]/(r*r*r);						

			Planet.position[0] += h*Planet.velocity[0]; 
			Planet.velocity[0] += h*a_x; 

			Planet.position[1] += h*Planet.velocity[1]; 
			Planet.velocity[1] += h*a_y; 

			ofile << Planet.position[0] << setw(20) << Planet.position[1] << endl;
		}
	
		//ostringstream os; 
		//os << "Output/EarthSun_"<< method  <<".txt"; 

		//string outfile = os.str(); 
		//ofile.open(outfile.c_str()); 

		// Write to file 
		ofile.close();

	}
	//print_position("Output_EarthSun_FE.txt",x,y); 


}  


void solver::VelocityVerlet(int integration_points, double time){

	// Constants & vectors 

	int n = integration_points; 	
	double h = time/((double) n);

	for(int i = 0; i<total_planets; i++){	 

		ofile.open("Output/EarthSun_VV.txt"); 
		planet &Planet = all_planets[i]; 

		// Initial conditions 

		double a_x; double a_y; double a_x_new; double a_y_new;  		 
		double r_0 = sqrt(Planet.position[0]*Planet.position[0]  + Planet.position[1]*Planet.position[1]); 
		a_x = - 4.0*M_PI*M_PI*Planet.position[0]/(r_0*r_0*r_0); 
		a_y = - 4.0*M_PI*M_PI*Planet.position[1]/(r_0*r_0*r_0); 

		// Solve with Velocity Verlet algo

		for(int i = 0; i<n-1; i++){

			Planet.position[0] += h*Planet.velocity[0] + h*h/2.0*a_x; 
			Planet.position[1] += h*Planet.velocity[1] + h*h/2.0*a_y; 

			r_0 = sqrt(Planet.position[0]*Planet.position[0]  + Planet.position[1]*Planet.position[1]); 
			a_x_new = - 4.0*M_PI*M_PI*Planet.position[0]/(r_0*r_0*r_0); 
			a_y_new = - 4.0*M_PI*M_PI*Planet.position[1]/(r_0*r_0*r_0); 
		
			Planet.velocity[0] += h/2.0*(a_x_new + a_x); 
			Planet.velocity[1] += h/2.0*(a_y_new + a_y); 

			a_x = a_x_new; a_y = a_y_new; 

			ofile << Planet.position[0] << setw(20) << Planet.position[1] << endl;

		}

		ofile.close(); 

	}
	
}

/*
void solver::print_position(string outname, vec x, vec y){
	ofstream ofile; 
	//ostringstream os; 
	//os << "Output/EarthSun_"<< method  <<".txt"; 

	//string outfile = os.str(); 
	//ofile.open(outfile.c_str()); 

	ofile.open(outname); 

	for(int i = 0; i < n; i++){
		ofile << x[i] << setw[20] << y[i] << endl; // << setw(20) << v_x(i) << setw(20) << v_y(i) << endl; 						
	}	

	ofile.close();
}
*/

 
