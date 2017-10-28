#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <string>
#include <iomanip>
#include <map>
#include "time.h"
#include "planet.h"
#include "solver.h"


using namespace std;

int main(int argc, char *argv[]){ 

	// Define some constants 

	string problem = argv[1]; // test, binary, three-body, solar system, mercury 
	int n = atoi(argv[2]);
	double t = atof(argv[3]);
	double variable; 
	if(problem == "Binary" or problem == "ThreeBody") variable = atof(argv[4]); 
	double yr = 365.25; // one year (in days)     
	double M_sun = 2.0E30; 

	// Test of the algorithms and escape velocity (c and d)  
	
	if(problem == "Binary"){
		planet test_earth(0.000003, 1.0, 0.0, 0.0, 2*M_PI, "Earth"); 	
		planet test_sun(1.0, 0, 0, 0, 0, "Sun"); 

		cout << "----------------------" << endl; 
		cout << "Initial values: " << endl; 
		cout << "Kinetic energy: " << test_earth.KineticEnergy() << endl; 
		cout << "Potential energy: " << test_earth.PotentialEnergy(test_sun) << endl; 
		cout << "Angular momentum: " << test_earth.AngularMomentum() << endl; 
		cout << "----------------------" << endl; 

		solver Test_FE;
		Test_FE.addPlanet(test_sun); 
		Test_FE.addPlanet(test_earth);

		clock_t FE_start, FE_finish; 
		FE_start = clock();  
		Test_FE.ForwardEuler(n, t, 1);  
		FE_finish = clock(); 

		solver Test_VV; 
		Test_VV.addPlanet(test_sun); 
		Test_VV.addPlanet(test_earth);	

		clock_t VV_start, VV_finish; 
		VV_start = clock(); 
		Test_VV.VelocityVerlet(n, t, 1); 
		VV_finish = clock(); 

		cout << "----------------------" << endl; 

		double FE_time = (FE_finish - FE_start)/((double) CLOCKS_PER_SEC); 
		double VV_time = (VV_finish - VV_start)/((double) CLOCKS_PER_SEC); 	

		cout << "Time with FE: " << FE_time << endl; 
		cout << "Time with VV: " << VV_time << endl; 
	}


	if(problem == "ThreeBody"){

		double v_e = 2*M_PI; double v_j = 2*M_PI*5.2/11.86; 
		double m_s = 1.0; double m_e = 0.000003; double m_j = variable*0.001; 

		planet sun(m_s, 0, 0, 0, 0, "Sun"); 
		planet earth(m_e, 1.0, 0.0, 0.0, v_e, "Earth"); 
		planet jupiter(m_j, 5.2, 0.0, 0.0,  v_j, "Jupiter"); 
			
		solver ThreeBody; 
		ThreeBody.addPlanet(sun);
		ThreeBody.addPlanet(earth); 
		ThreeBody.addPlanet(jupiter);  
		//ThreeBody.VelocityVerlet(n, t, 1);

		
		// Calculate center-of-mass (only necessary for x-dir with the given initial conditions) 
		double M_tot = m_s + m_e + m_j; 
		double cm_x = 1/M_tot*(m_e*1.0 + m_j*5.2); 
		//cout << cm_x << endl; 
		//double cm_y = 1/M_tot*(0.0 + 0.0); 		

		// Updating 

		sun.position[0] = sun.position[0] - cm_x; 
		earth.position[0] = earth.position[0] - cm_x; 
		jupiter.position[0] = jupiter.position[0] - cm_x; 

		sun.velocity[1] = -(earth.mass*earth.velocity[1] + jupiter.mass*jupiter.velocity[1]);		

		solver ThreeBody_CM; 
		ThreeBody_CM.addPlanet(sun); 
		ThreeBody_CM.addPlanet(earth); 
		ThreeBody_CM.addPlanet(jupiter);
		ThreeBody_CM.VelocityVerlet(n,t,1); 

 
	}

	if(problem == "SolarSystem"){

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

	// Intialize planets 

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
	SolarSystem.VelocityVerlet(n, t, 1); 

}

if(problem == "Mercury"){

	// Perihelion precession (sun and mercury)

	planet sun(1.0, 0, 0, 0, 0, "Sun"); 
	planet mercury(1.65E-7, 0.3075, 0.0, 0.0, 12.44, "Mercury"); 	

	solver Mercury_GR(1); 
	Mercury_GR.addPlanet(sun); 
	Mercury_GR.addPlanet(mercury); 
	Mercury_GR.VelocityVerlet(n, t, 0); 	

	solver Mercury_Newton; 
	Mercury_Newton.addPlanet(sun); 
	Mercury_Newton.addPlanet(mercury); 
	Mercury_Newton.VelocityVerlet(n, t, 0); 	
	

}

}
