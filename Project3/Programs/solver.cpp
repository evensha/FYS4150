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
	beta = 2.0; 
	RelCorr = 0; 

}


solver::solver(double b){

	total_planets = 0; 
	mass = 0; 
	radius = 0; 
	G = 4*M_PI*M_PI; 
	beta = b; 
	RelCorr = 0; 

}

solver::solver(int withGR){

	total_planets = 0; 
	mass = 0; 
	radius = 0; 
	G = 4*M_PI*M_PI; 
	beta = 2.0; 
	RelCorr = withGR; 

}


/*
solver::solver(double r){

	total_planets = 0; 
	mass = 0; 
	radius = r; 
	G = 4*M_PI*M_PI; 
	beta = 2.0; 

}
*/
	

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

	double F_x=0; double F_y=0; double F_z=0;  	
	double a_x=0; double a_y=0; double a_z=0; 

	for(int i = 0; i<n; i++){

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

			// Update positions and velocities 	

			Planet.position[0] += h*Planet.velocity[0]; 
			Planet.velocity[0] += h*a_x; 

			Planet.position[1] += h*Planet.velocity[1]; 
			Planet.velocity[1] += h*a_y; 

			Planet.position[2] += h*Planet.velocity[2]; 
			Planet.velocity[2] += h*a_z; 

			F_x = 0; F_y = 0; F_z = 0; // reset forces 

		}

		PrintPositions(); 		

	}

	ofile.close();

}  


void solver::VelocityVerlet(int integration_points, double final_time){

	// Constants & vectors 

	int n = integration_points; 	
	double h = final_time/((double) n);	

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

	double F_x=0; double F_y=0; double F_z=0; double F_x_new=0; double F_y_new=0; double F_z_new=0;  	
	double a_x=0; double a_y=0; double a_z=0; double a_x_new=0; double a_y_new=0; double a_z_new=0;  		
	double smallest_r; double x_peri; double y_peri; int smallest_r_time;  

	for(int i = 1; i<n; i++){    // Start loop over time steps
	
		for(int j = 0; j<total_planets; j++){	 
			planet &Planet = all_planets[j]; 
			if(i==1) cout << Planet.name << endl; 
			if(Planet.name == "Sun") continue;  // sun fixed -> don't bother with the calculations

			for(int k = 0; k<total_planets; k++){  // calculate gravitational forces on the planet 
				if(k != j){
					planet &other = all_planets[k]; 
					GravitationalForce(Planet, other, F_x, F_y, F_z, beta, RelCorr); 
				}
			}

			double m = Planet.mass;  
			a_x = F_x/m; a_y = F_y/m; a_z=F_z/m;  // acceleration 

			double r = Planet.Distance(all_planets[0]); 

			// Solve with Velocity Verlet algo

			Planet.position[0] += h*Planet.velocity[0] + 0.5*h*h*a_x;  // new positions 
			Planet.position[1] += h*Planet.velocity[1] + 0.5*h*h*a_y; 
			Planet.position[2] += h*Planet.velocity[2] + 0.5*h*h*a_z; 

			if( Planet.Distance(all_planets[0]) <= r ){ 
				smallest_r = Planet.Distance(all_planets[0]); 
				smallest_r_time = i;
				x_peri = Planet.position[0]; 
				y_peri = Planet.position[1]; 
			} 
			/*
			if(Planet.Distance(all_planets[0]) > 0.3074999 && Planet.Distance(all_planets[0]) < 0.3075001 && i>0.96*n){ 
				cout << "Approximately at perihelion" << endl; 
				cout << Planet.position[0] << endl; 
				cout << Planet.position[1] << endl; 
				cout << "tan t_p = " << Planet.position[1]/Planet.position[0] << endl; 
				cout << i << endl; 
				cout << "------------------" << endl; 
			} 
			*/
			for(int k = 0; k<total_planets; k++){
				if( k != j ){
					planet &other = all_planets[k]; 
					GravitationalForce(Planet, other, F_x_new, F_y_new, F_z_new, beta, RelCorr); 	// new forces
				}
			}	

			a_x_new = F_x_new/m; a_y_new = F_y_new/m; a_z_new = F_z_new/m;  		
		
			Planet.velocity[0] += 0.5*h*(a_x_new + a_x);  // velocities 
			Planet.velocity[1] += 0.5*h*(a_y_new + a_y); 
			Planet.velocity[2] += 0.5*h*(a_z_new + a_z); 

			a_x = a_x_new; a_y = a_y_new; a_z = a_z_new;   // update acceleration 
			F_x = 0; F_y = 0; F_z = 0; F_x_new = 0; F_y_new = 0; F_z_new = 0;   // reset forces 	

		}

		//PrintPositions();  // print updated positions
		//time += h; 

	}

	cout << "Smallest r and time: " << endl; 
	cout << smallest_r << endl; 
	cout << x_peri << endl; 
	cout << y_peri << endl; 
	cout << "tan t_p = " << y_peri/x_peri << endl; 
	cout << smallest_r_time << endl; 
	ofile.close(); 
	
}


void solver::GravitationalForce(planet &Planet, planet &other, double &F_x, double &F_y, double &F_z, double beta, int RelCorr){

	double dist_x =  other.position[0] - Planet.position[0]; 
	double dist_y =  other.position[1] - Planet.position[1];  
	double dist_z =  other.position[2] - Planet.position[2]; 
	double r = Planet.Distance(other); 
	//cout << r << endl; 
	double m = Planet.mass; double M = other.mass; 

	double x = Planet.position[0]; double y = Planet.position[1]; double z = Planet.position[2]; 
	double v_x = Planet.velocity[0]; double v_y = Planet.velocity[1]; double v_z = Planet.velocity[2]; 

	double L_x = y*v_z - z*v_y; 
	double L_y = -(x*v_z - z*v_x);
	if( L_x != 0 or L_y != 0) cout << "Something wrong!" << endl;  
	double L_z = x*v_y - y*v_x; 
	double l = sqrt(L_x*L_x + L_y*L_y + L_z*L_z); 

	//double l = Planet.AngularMomentum(); 
	double c = 63284.9; // speed of light in AU/year 

	double GR = 0.0; 
	if(RelCorr == 1){
		GR = 3.0*l*l/(r*r*c*c); 
	}

	//cout << "---------------" << endl; 
	//cout << "GR correction: " << GR << endl; 

	F_x += G*m*M*dist_x/(r*r*r)*(1.0 + GR); 
	F_y += G*m*M*dist_y/(r*r*r)*(1.0 + GR); 
	F_z += G*m*M*dist_z/(r*r*r)*(1.0 + GR); 		  

	//cout << "Force in x-dir: " << G*m*M*dist_x/(r*pow(r,beta)) << endl; 
	//cout << "Corrected force: " << F_x << endl; 

}


void solver::PrintPositions(){
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
	for(int i = 0; i<total_planets; i++){
		if( i == total_planets-1) ofile << all_planets[i].name; 
		else ofile << all_planets[i].name << setw(60);
	}
	ofile << endl; 
}



 
