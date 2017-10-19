#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <string>
#include <iomanip>
#include <map>
//#include <armadillo>
#include "time.h"
#include "planet.h"
#include "solver.h"


using namespace std;

int main(int argc, char *argv[]){ 

	// Define some constants 

	int n = atoi(argv[1]);
	int withGR = atoi(argv[3]); 
	double t = atof(argv[2]); 
	double yr = 365.25; // one year (in days)     
	double M_sun = 2.0E30; 

	// Make maps for planet data 

	map<int, string> names;
	map<string, double> mass; 
	map<string, double> x;
	map<string, double> y; 
	map<string, double> z; 
	map<string, double> v_x;
	map<string, double> v_y;
	map<string, double> v_z;  

	// Extract planet data from file and store in the maps

	ifstream infile("Planet_data.txt"); 

	string Planet; double M; double X; double Y; double Z; double V_X; double V_Y; double V_Z; 

	int i = 0; 
	while (infile >> Planet >> M >> X >> Y >> Z >> V_X >> V_Y >> V_Z  ){
		names[i] = Planet; 
		mass[Planet] = M/M_sun; x[Planet] = X; y[Planet] = Y; z[Planet] = Z; v_x[Planet] = V_X*yr; v_y[Planet] = V_Y*yr; v_z[Planet] = V_Z*yr;  
		i += 1; 
	} 
 
	infile.close(); 

	// Intialize planets (for now only in 2D) 

	planet sun(mass["Sun"], x["Sun"], y["Sun"], z["Sun"], v_x["Sun"], v_y["Sun"], v_z["Sun"], "Sun"); 
	planet mercury(mass["Mercury"], x["Mercury"], y["Mercury"], z["Mercury"], v_x["Mercury"], v_y["Mercury"], v_z["Mercury"], "Mercury");
	planet venus(mass["Venus"], x["Venus"], y["Venus"], z["Venus"], v_x["Venus"], v_y["Venus"], v_z["Venus"], "Venus"); 
	planet earth(mass["Earth"], x["Earth"], y["Earth"], z["Earth"], v_x["Earth"], v_y["Earth"], v_z["Earth"], "Earth"); 
	planet mars(mass["Mars"], x["Mars"], y["Mars"], z["Mars"], v_x["Mars"], v_y["Mars"], v_z["Mars"], "Mars");  
	planet jupiter(mass["Jupiter"], x["Jupiter"], y["Jupiter"], z["Jupiter"], v_x["Jupiter"], v_y["Jupiter"], v_z["Jupiter"], "Jupiter"); 
	planet saturn(mass["Saturn"], x["Saturn"], y["Saturn"], z["Saturn"], v_x["Saturn"], v_y["Saturn"], v_z["Saturn"], "Saturn"); 
	planet uranus(mass["Uranus"], x["Uranus"], y["Uranus"], z["Uranus"], v_x["Uranus"], v_y["Uranus"], v_z["Uranus"], "Uranus"); 
	planet neptun(mass["Neptun"], x["Neptun"], y["Neptun"], z["Neptun"], v_x["Neptun"], v_y["Neptun"], v_z["Neptun"], "Neptun"); 
	planet pluto(mass["Pluto"], x["Pluto"], y["Pluto"], z["Pluto"], v_x["Pluto"], v_y["Pluto"], v_z["Pluto"], "Pluto"); 
	//planet earth(0.000003, 1.0, 0.0, 0.0, v_start*2*M_PI, "Earth"); 	

	// Print some initial properties of earth

	cout << "----------------------" << endl; 
	cout << "Initial values: " << endl; 
	cout << "Kinetic energy: " << earth.KineticEnergy() << endl; 
	cout << "Potential energy: " << earth.PotentialEnergy(4*M_PI*M_PI) << endl; 
	cout << "Angular momentum: " << earth.AngularMomentum() << endl; 
	cout << "----------------------" << endl; 

	// Initialize solver class, add planets and solve the problem 

	// Binary system (earth and sun) 
/*
	solver Binary_FE;
	Binary_FE.addPlanet(sun); 
	Binary_FE.addPlanet(earth);

	clock_t FE_start, FE_finish; 
	FE_start = clock();  
	Binary_FE.ForwardEuler(n, t);  
	FE_finish = clock(); 

	solver Binary_VV; 
	Binary_VV.addPlanet(sun); 
	Binary_VV.addPlanet(earth);	

	clock_t VV_start, VV_finish; 
	VV_start = clock(); 
	Binary_VV.VelocityVerlet(n, t); 
	VV_finish = clock(); 

	// Three-body system (earth, sun and jupiter) 

	sun.velocity[0] = -( v_x["Earth"]*mass["Earth"] + v_x["Jupiter"]*mass["Jupiter"] );  
	sun.velocity[1] = -( v_y["Earth"]*mass["Earth"] + v_y["Jupiter"]*mass["Jupiter"] );  

	solver ThreeBody; 
	ThreeBody.addPlanet(sun);
	ThreeBody.addPlanet(earth); 
	ThreeBody.addPlanet(jupiter);  
	ThreeBody.VelocityVerlet(n, t); 

	// Solar system (all planets) 

	sun.velocity[0] = 0.0; sun.velocity[1] = 0.0; string planet;  // reset position and velocity of sun  
	sun.position[0] = x["Sun"]; sun.position[1] = y["Sun"]; 
	for(int i = 1; i<10; i++){
		planet = names[i]; 
		sun.velocity[0] -= v_x[planet]*mass[planet]; 
		sun.velocity[1] -= v_y[planet]*mass[planet]; 
	}

	solver SolarSystem; 
	SolarSystem.addPlanet(sun); 	
	SolarSystem.addPlanet(mercury); 	
	SolarSystem.addPlanet(venus);
	SolarSystem.addPlanet(earth);	
	SolarSystem.addPlanet(mars);
	SolarSystem.addPlanet(jupiter);
	SolarSystem.addPlanet(saturn);
	SolarSystem.addPlanet(uranus);
	SolarSystem.addPlanet(neptun);
	SolarSystem.addPlanet(pluto);
	SolarSystem.VelocityVerlet(n, t); 
*/
	// Perihelion precession (sun and mercury)

	mercury.position[0] = 0.3075; mercury.position[1] = 0.0; mercury.position[2] = 0.0; 
	mercury.velocity[0] = 0.0; mercury.velocity[1] = 12.44; mercury.velocity[2] = 0.0; 
	 // reset position and velocity of sun  
	sun.position[0] = 0.0; sun.position[1] = 0.0; sun.position[2] = 0.0; 
	sun.velocity[0] = 0.0; sun.velocity[1] = 0.0; sun.velocity[2] = 0.0; 

	solver Mercury(withGR); 
	Mercury.addPlanet(sun); 
	Mercury.addPlanet(mercury); 
	Mercury.VelocityVerlet(n, t); 	
/*
	cout << "----------------------" << endl; 

	double FE_time = (FE_finish - FE_start)/((double) CLOCKS_PER_SEC); 
	double VV_time = (VV_finish - VV_start)/((double) CLOCKS_PER_SEC); 	

	cout << "Time with FE: " << FE_time << endl; 
	cout << "Time with VV: " << VV_time << endl; 
*/
}
