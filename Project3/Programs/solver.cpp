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
	// Intialize solver with no specifications 

	total_planets = 0; 
	mass = 0; 
	G = 4*M_PI*M_PI; 
	beta = 2.0; 
	RelCorr = 0; 

}


solver::solver(double b){
	// Initialize solver with a double; used to alter the form of the gravitational force 

	total_planets = 0; 
	mass = 0; 
	G = 4*M_PI*M_PI; 
	beta = b; 
	RelCorr = 0; 

}

solver::solver(int withGR){
	// Initialize solver with an int; used to indicate whether or not to calculate relativistic corrections 

	total_planets = 0; 
	mass = 0; 
	G = 4*M_PI*M_PI; 
	beta = 2.0; 
	RelCorr = withGR; 

}
	

void solver::addPlanet(planet newplanet){ 
	// Add a planet to the solver class 

	total_planets += 1; 
	all_planets.push_back(newplanet); 
	mass += newplanet.mass; 

} 


void solver::ForwardEuler(int integration_points, double time, int withOutput){
	// Solve the equations of motion for the added planets using the forward Euler algorithm for given time period and number of 
	// integration points  

	int n = integration_points; 	
	double h = time/((double) n);

	// Make outfile 

	string problem; 
	if(withOutput == 1){
		ostringstream os; 
		if(total_planets == 2) problem = "Binary"; 
		else if(total_planets == 3) problem = "ThreeBody"; 
		else if(total_planets == 10) problem = "SolarSystem"; 
		else problem = "Planets"; 

		if(total_planets == 2 && all_planets[1].name == "Mercury") problem = "Mercury_perihelion"; 

		if(RelCorr == 1) os << "Output/" << problem << "_FE_withGR" << ".txt";
		else if( beta != 2 ) os << "Output/" << problem << "_FE_beta=" << beta << ".txt";
		else os << "Output/" << problem << "_FE" << ".txt";  

		string outfile = os.str(); 
		ofile.open(outfile.c_str()); 

		PrintNames(); 
		PrintPositions();  // print initial positions 
	}


	// Initialize forces and acceleration

	double F_x=0; double F_y=0; double F_z=0;  	
	double a_x=0; double a_y=0; double a_z=0; 

 	// Loop over time

	for(int i = 0; i<n; i++){  

		// Loop over planets

		for(int j = 0; j<total_planets; j++){	  

			planet &Planet = all_planets[j]; 
			if(Planet.name == "Sun") continue; 

			for(int k = 0; k<total_planets; k++){  // calculate gravitational forces on the planet 
				if(k != j){
					planet &other = all_planets[k]; 
					GravitationalForce(Planet, other, F_x, F_y, F_z, beta, RelCorr); 
				}
			}

			double m = Planet.mass;  
			a_x = F_x/m; a_y = F_y/m; a_z = F_z/m; // acceleration 	

			// Update positions and velocities according to the forward Euler method 

			Planet.position[0] += h*Planet.velocity[0]; 
			Planet.velocity[0] += h*a_x; 

			Planet.position[1] += h*Planet.velocity[1]; 
			Planet.velocity[1] += h*a_y; 

			Planet.position[2] += h*Planet.velocity[2]; 
			Planet.velocity[2] += h*a_z; 

			F_x = 0; F_y = 0; F_z = 0; // Reset forces before jumping to the next planet  

		}

		PrintPositions(); // Write positions to file after each time step   		

	}

	ofile.close();

	if(problem == "Binary"){
		cout << "----------------------" << endl; 
		cout << "After forward Euler:" << endl; 
		cout << "Kinetic energy: " << all_planets[1].KineticEnergy() << endl; 
		cout << "Potential energy: " << all_planets[1].PotentialEnergy(all_planets[0]) << endl; 
		cout << "Angular momentum: " << all_planets[1].AngularMomentum() << endl; 
		cout << "----------------------" << endl; 
	}

}  


void solver::VelocityVerlet(int integration_points, double final_time, int withOutput){
	// Solve the equations of motion for the added planets using the velocity Verlet algorithm for given time period and number of 
	// integration points. Also indicate whether or not you want to produce an output file.   

	int n = integration_points; 	
	double h = final_time/((double) n);	

	double h_half = 0.5*h; 
	double h2_half = 0.5*h*h; 

	// Make outfile (if we want output) 

	string problem; 
	if(withOutput == 1){
		ostringstream os; 
		if(total_planets == 2) problem = "Binary"; 
		else if(total_planets == 3) problem = "ThreeBody"; 
		else if(total_planets == 10) problem = "SolarSystem"; 
		else problem = "Planets"; 

		if(total_planets == 2 && all_planets[1].name == "Mercury") problem = "Mercury_perihelion"; 

		if(RelCorr == 1) os << "Output/" << problem << "_VV_withGR" << ".txt";
		else if( beta != 2 ) os << "Output/" << problem << "_VV_beta=" << beta << ".txt";
		else os << "Output/" << problem << "_VV" << ".txt";  

		string outfile = os.str(); 
		ofile.open(outfile.c_str()); 

		PrintNames(); 
		PrintPositions();  // print initial positions 
	}

	// Initialize forces and acceleration

	double F_x=0; double F_y=0; double F_z=0; double F_x_new=0; double F_y_new=0; double F_z_new=0;  	
	double a_x=0; double a_y=0; double a_z=0; double a_x_new=0; double a_y_new=0; double a_z_new=0;  	
	
	double smallest_r; double x_peri; double y_peri; int smallest_r_time;  // quantities related to perihelion   

	// Loop over time   

	for(int i = 1; i<n; i++){   
	
		// Loop over planets 

		for(int j = 0; j<total_planets; j++){	 
			planet &Planet = all_planets[j]; 
			//if(Planet.name == "Sun") continue;  // Sun fixed -> don't bother with the calculations

			for(int k = 0; k<total_planets; k++){  // Calculate gravitational forces on the planet 
				if(k != j){
					planet &other = all_planets[k]; 
					GravitationalForce(Planet, other, F_x, F_y, F_z, beta, RelCorr); 
				}
			}

			double m = Planet.mass;  
			a_x = F_x/m; a_y = F_y/m; a_z=F_z/m;  // Acceleration 

			double r = Planet.Distance(all_planets[0]); 

			// Solve with Velocity Verlet algo

			Planet.position[0] += h*Planet.velocity[0] + h2_half*a_x;  // new positions 
			Planet.position[1] += h*Planet.velocity[1] + h2_half*a_y; 
			Planet.position[2] += h*Planet.velocity[2] + h2_half*a_z; 

			if(withOutput == 0){  // study perihelion precession for the binary system 
				if( Planet.Distance(all_planets[0]) <= r ){   // update perihelion position 
					smallest_r = Planet.Distance(all_planets[0]); 
					smallest_r_time = i;
					x_peri = Planet.position[0]; 
					y_peri = Planet.position[1]; 
				} 
			}
	
			for(int k = 0; k<total_planets; k++){
				if( k != j ){
					planet &other = all_planets[k]; 
					GravitationalForce(Planet, other, F_x_new, F_y_new, F_z_new, beta, RelCorr); 	// New forces
				}
			}	

			a_x_new = F_x_new/m; a_y_new = F_y_new/m; a_z_new = F_z_new/m;  		
		
			Planet.velocity[0] += h_half*(a_x_new + a_x);  // New velocities 
			Planet.velocity[1] += h_half*(a_y_new + a_y); 
			Planet.velocity[2] += h_half*(a_z_new + a_z); 

			a_x = a_x_new = a_y = a_y_new = a_z = a_z_new = 0; 
			//a_x = a_x_new; a_y = a_y_new; a_z = a_z_new;   // Update acceleration 
			F_x = F_y = F_z = F_x_new = F_y_new = F_z_new = 0;   // Reset forces 	

		}

		if(withOutput == 1) PrintPositions();  // Print updated positions 

	}

	if(withOutput == 0){  // Print out perhelion info for the binary system 
		cout << "Smallest r and time: " << endl; 
		cout << smallest_r << endl; 
		cout << x_peri << endl; 
		cout << y_peri << endl; 
		cout << "tan t_p = " << y_peri/x_peri << endl; 
		cout << smallest_r_time << endl; 
	}

	ofile.close(); 

	if(problem == "Binary"){
		cout << "----------------------" << endl; 
		cout << "After velocity Verlet:" << endl; 
		cout << "Kinetic energy: " << all_planets[1].KineticEnergy() << endl; 
		cout << "Potential energy: " << all_planets[1].PotentialEnergy(all_planets[0]) << endl; 
		cout << "Angular momentum: " << all_planets[1].AngularMomentum() << endl; 
		cout << "----------------------" << endl; 
	}
	
}


void solver::GravitationalForce(planet &Planet, planet &other, double &F_x, double &F_y, double &F_z, double beta, int RelCorr){
	// Calculate the gravitational force between two planets 

	double dist_x =  other.position[0] - Planet.position[0]; 
	double dist_y =  other.position[1] - Planet.position[1];  
	double dist_z =  other.position[2] - Planet.position[2]; 
	double r = Planet.Distance(other); 
	double m = Planet.mass; double M = other.mass; 
	
	double l = Planet.AngularMomentum(); 
	double c = 63284.9; // Speed of light in AU/year 

	double GR = 0.0; 
	if(RelCorr == 1){
		GR = 3.0*l*l/(r*r*c*c); 
	}

	F_x += G*m*M*dist_x/(pow(r,beta)*r)*(1.0 + GR); 
	F_y += G*m*M*dist_y/(pow(r,beta)*r)*(1.0 + GR); 
	F_z += G*m*M*dist_z/(pow(r,beta)*r)*(1.0 + GR); 		  

}



void solver::PrintPositions(){
	// Print positions of all your planets to the outfile 

	for(int i = 0; i<total_planets; i++){
		if( i == total_planets-1){
			ofile << setprecision(8) << all_planets[i].position[0] << setw(20) 
			<< all_planets[i].position[1] << setw(20) << all_planets[i].position[2];
		}
		else{
			ofile << setprecision(8) << all_planets[i].position[0] << setw(20) << 
			all_planets[i].position[1] << setw(20) << all_planets[i].position[2] << setw(20); 
		}
	}
	ofile << endl; 

}


void solver::PrintNames(){
	// Print the names of your planets to the outfile 

	for(int i = 0; i<total_planets; i++){
		if( i == total_planets-1) ofile << all_planets[i].name; 
		else ofile << all_planets[i].name << setw(60);
	}
	ofile << endl; 
}


 
