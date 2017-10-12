#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <string>
#include <iomanip>
//#include <armadillo>
#include "time.h"
#include "planet.h"
#include "solver.h"


using namespace std;

int main(int argc, char *argv[]){ 

	// Define some constants 

	int n = atoi(argv[1]);
	double t = atof(argv[2]); 

	double v_start = atof(argv[3]);   

	double year = 365.25; 
	double M_sun = 1.0; 
	double M_earth = 0.000003;
	double M_mars = 3.3E-7;  
	double M_jupiter = 0.00095;   

	// Intialize planets 	

	planet sun(M_sun, 0, 0, 0, 0, "Sun"); 
	planet earth(M_earth, 9.734E-01, 2.4198E-01, -4.349E-03*year, 1.665E-02*year, "Earth"); 
	planet mars(M_mars, -1.5158E+00, 6.9015E-01, -5.2319E-03*year, -1.1556E-02*year, "Mars"); 
	planet jupiter(M_jupiter, -4.623E+00, -2.860E+00, 3.881E-03*year, -6.058E-03*year, "Jupiter"); 
	//planet earth(0.000003, 1.0, 0.0, 0.0, v_start*2*M_PI, "Earth"); 	

	cout << "----------------------" << endl; 
	cout << "Initial values: " << endl; 
	cout << "Kinetic energy: " << earth.KineticEnergy() << endl; 
	cout << "Potential energy: " << earth.PotentialEnergy(4*M_PI*M_PI) << endl; 
	cout << "Angular momentum: " << earth.AngularMomentum() << endl; 

	// Initialize solver class, add planets and solve the problem 

	clock_t FE_start, FE_finish; 
	FE_start = clock(); 

	solver EarthSun_FE(1.0); 
	EarthSun_FE.addPlanet(earth); 
	EarthSun_FE.ForwardEuler(n, t);  

	FE_finish = clock(); 

	clock_t VV_start, VV_finish; 
	VV_start = clock(); 

	solver Planets_VV(1.0); 
	Planets_VV.addPlanet(sun);
	Planets_VV.addPlanet(earth); 
	Planets_VV.addPlanet(mars); 
	Planets_VV.addPlanet(jupiter);  
	Planets_VV.VelocityVerlet(n, t); 

	VV_finish = clock(); 

	solver EarthSun_VV(1.0); 
	EarthSun_VV.addPlanet(sun); 
	EarthSun_VV.addPlanet(earth); 
	EarthSun_VV.VelocityVerlet(n,t); 

	cout << "----------------------" << endl; 

	double FE_time = (FE_finish - FE_start)/((double) CLOCKS_PER_SEC); 
	double VV_time = (VV_finish - VV_start)/((double) CLOCKS_PER_SEC); 	

	cout << "Time with FE: " << FE_time << endl; 
	cout << "Time with VV: " << VV_time << endl; 

}
