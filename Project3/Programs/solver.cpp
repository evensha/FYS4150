#include "solver.h"
#include "planet.h"
#include <iostream> 
#include <cmath>
#include "time.h"
#include <cstdlib> 
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

using namespace std; 

ofstream ofile; 

solver::solver(){

	total_planets = 0; 
	mass = 0; 
	radius = 0; 
	G = 4*M_PI*M_PI; 

}

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

	// Make outfile 

	ostringstream os; 
	string problem; 
	if(total_planets == 2) problem = "Binary"; 
	else if(total_planets == 3) problem = "ThreeBody"; 
	else if(total_planets == 10) problem = "SolarSystem"; 
	else problem = "Planets"; 
	os << "Output/" << problem << "_FE" << ".txt"; 

	string outfile = os.str(); 
	ofile.open(outfile.c_str()); 

	PrintNames(); 
	PrintPositions();  // print initial positions 

	// Initialize forces and acceleration

	double F_x=0; double F_y=0;   	
	double a_x=0; double a_y=0;

	for(int i = 0; i<n; i++){

		for(int j = 0; j<total_planets; j++){	 

			planet &Planet = all_planets[j]; 
			if(Planet.name == "Sun") continue; 

			for(int k = 0; k<total_planets; k++){  // calculate gravitational forces on the planet 
				if(k != j){
					planet &other = all_planets[k]; 
					GravitationalForce(Planet, other, F_x, F_y); 
				}
			}

			double m = Planet.mass;  
			a_x = F_x/m; a_y = F_y/m;  // acceleration 	

			// Update positions and velocities 	

			Planet.position[0] += h*Planet.velocity[0]; 
			Planet.velocity[0] += h*a_x; 

			Planet.position[1] += h*Planet.velocity[1]; 
			Planet.velocity[1] += h*a_y; 

			F_x = 0; F_y = 0; // reset forces 

		}

		PrintPositions(); 		

	}

	ofile.close();

}  


void solver::VelocityVerlet(int integration_points, double time){

	// Constants & vectors 

	int n = integration_points; 	
	double h = time/((double) n);	

	// Make outfile 

	ostringstream os; 
	string problem; 
	if(total_planets == 2) problem = "Binary"; 
	else if(total_planets == 3) problem = "ThreeBody"; 
	else if(total_planets == 10) problem = "SolarSystem"; 
	else problem = "Planets"; 
	os << "Output/" << problem << "_VV" << ".txt";  

	string outfile = os.str(); 
	ofile.open(outfile.c_str()); 

	PrintNames(); 
	PrintPositions();  // print initial positions 

	// Initialize forces and acceleration

	double F_x=0; double F_y=0; double F_x_new=0; double F_y_new=0;  	
	double a_x=0; double a_y=0; double a_x_new=0; double a_y_new=0;  		
 
	for(int i = 0; i<n; i++){    // Start loop over time steps
	
		for(int j = 0; j<total_planets; j++){	 
			planet &Planet = all_planets[j]; 
			if(Planet.name == "Sun") continue;  // sun fixed -> don't bother with the calculations

			for(int k = 0; k<total_planets; k++){  // calculate gravitational forces on the planet 
				if(k != j){
					planet &other = all_planets[k]; 
					GravitationalForce(Planet, other, F_x, F_y); 
				}
			}

			double m = Planet.mass;  
			a_x = F_x/m; a_y = F_y/m;  // acceleration 

			// Solve with Velocity Verlet algo

			Planet.position[0] += h*Planet.velocity[0] + h*h/2.0*a_x;  // new positions 
			Planet.position[1] += h*Planet.velocity[1] + h*h/2.0*a_y; 

			for(int k = 0; k<total_planets; k++){
				if( k != j ){
					planet &other = all_planets[k]; 
					GravitationalForce(Planet, other, F_x_new, F_y_new); 	// new forces
				}
			}	
			a_x_new = F_x_new/m; a_y_new = F_y_new/m; 		
		
			Planet.velocity[0] += h/2.0*(a_x_new + a_x);  // velocities 
			Planet.velocity[1] += h/2.0*(a_y_new + a_y); 

			a_x = a_x_new; a_y = a_y_new;   // update acceleration 
			F_x = 0; F_y = 0; F_x_new = 0; F_y_new = 0;  // reset forces 	

		}

		PrintPositions();  // print updated positions

	}

	ofile.close(); 
	
}


void solver::GravitationalForce(planet &Planet, planet &other, double &F_x, double &F_y){

	double dist_x =  Planet.position[0] - other.position[0]; 
	double dist_y =  Planet.position[1] - other.position[1];  
	double r = Planet.Distance(other); 
	double m = Planet.mass; double M = other.mass; 

	F_x -= G*m*M*dist_x/(r*r*r); 
	F_y -= G*m*M*dist_y/(r*r*r);   

}


void solver::PrintPositions(){
	for(int i = 0; i<total_planets; i++){
		if( i == total_planets-1)	ofile << setprecision(5) << all_planets[i].position[0] << setw(20) << all_planets[i].position[1];
		else ofile << setprecision(5) << all_planets[i].position[0] << setw(20) << all_planets[i].position[1] << setw(20); 
	}
	ofile << endl; 

}


void solver::PrintNames(){
	for(int i = 0; i<total_planets; i++){
		if( i == total_planets-1) ofile << all_planets[i].name; 
		else ofile << all_planets[i].name << setw(40);
	}
	ofile << endl; 
}



 
