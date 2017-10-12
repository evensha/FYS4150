#include "solver.h"
#include "planet.h"
#include <iostream> 
#include <cmath>
#include "time.h"
#include <cstdlib> 
#include <fstream>
#include <iomanip>
#include <string>

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

		double a_x; double a_y; double r; 

		// Solve with Forward-Euler algo 

		for(int i = 0; i<n-1; i++){
			r = sqrt(Planet.position[0]*Planet.position[0] + Planet.position[1]*Planet.position[1]); 
			a_x = - 4.0*M_PI*M_PI*Planet.position[0]/(r*r*r); 
			a_y = - 4.0*M_PI*M_PI*Planet.position[1]/(r*r*r);						

			Planet.position[0] += h*Planet.velocity[0]; 
			Planet.velocity[0] += h*a_x; 

			Planet.position[1] += h*Planet.velocity[1]; 
			Planet.velocity[1] += h*a_y; 

			ofile << Planet.position[0] << setw(20) << Planet.position[1] << endl;
		}

		cout << "----------------------" << endl; 
		cout << "After Euler: " << endl; 
		cout << "Kinetic energy: " << Planet.KineticEnergy() << endl; 
		cout << "Potential energy: " << Planet.PotentialEnergy(G) << endl; 
		cout << "Angular momentum: " << Planet.AngularMomentum() << endl; 

		//ostringstream os; 
		//os << "Output/EarthSun_"<< method  <<".txt"; 

		//string outfile = os.str(); 
		//ofile.open(outfile.c_str()); 

		ofile.close();

	}
	//print_position("Output_EarthSun_FE.txt",x,y); 

}  


void solver::VelocityVerlet(int integration_points, double time){

	// Constants & vectors 

	int n = integration_points; 	
	double h = time/((double) n);

	for(int i = 0; i<total_planets-1; i++){	 

		ofile.open("Output/EarthSun_VV.txt"); 
		//ostringstream os; 
		//os << "Output/EarthSun_VV.txt"; 

		//string outfile = os.str(); 
		//ofile.open(outfile.c_str()); 

		//ofile.open(outname); 

		planet &Planet = all_planets[i]; 

		// Initial conditions 

		double F_x=0; double F_y=0; double F_x_new=0; double F_y_new=0;  	
		double a_x=0; double a_y=0; double a_x_new=0; double a_y_new=0;  		 

		GravitationalForce(Planet, all_planets[1], F_x, F_y); 
		double m = Planet.mass;  
		a_x = F_x/m; a_y = F_y/m; 

		// Solve with Velocity Verlet algo

		for(int i = 0; i<n-1; i++){

			Planet.position[0] += h*Planet.velocity[0] + h*h/2.0*a_x; 
			Planet.position[1] += h*Planet.velocity[1] + h*h/2.0*a_y; 

			GravitationalForce(Planet, all_planets[1], F_x_new, F_y_new); 	
			a_x_new = F_x_new/m; a_y_new = F_y_new/m; 		
		
			Planet.velocity[0] += h/2.0*(a_x_new + a_x); 
			Planet.velocity[1] += h/2.0*(a_y_new + a_y); 

			a_x = a_x_new; a_y = a_y_new; 
			F_x = 0; F_y = 0; F_x_new = 0; F_y_new = 0; 	

			ofile << Planet.position[0] << setw(20) << Planet.position[1] << endl;
	

		}

		cout << "----------------------" << endl; 
		cout << "After Verlet: " << endl; 
		cout << "Kinetic energy: " << Planet.KineticEnergy() << endl; 
		cout << "Potential energy: " << Planet.PotentialEnergy(G) << endl; 
		cout << "Angular momentum: " << Planet.AngularMomentum() << endl; 

		ofile.close(); 

	}
	
}


void solver::GravitationalForce(planet &Planet, planet &other, double &F_x, double &F_y){

	double dist_x = Planet.position[0]; double dist_y=Planet.position[1];  
	double r = Planet.Distance(other); 
	double m = Planet.mass; double M = other.mass; 

	F_x -= G*m*M*dist_x/(r*r*r); 
	F_y -= G*m*M*dist_y/(r*r*r);   

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

 
