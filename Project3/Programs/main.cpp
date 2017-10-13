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

	//double v_start = atof(argv[3]);   

	double year = 365.25; 
	double M_sun = 2.0E30;   
	double M_mercury = 3.3E23/M_sun; // masses relative to the sun mass
	double M_venus = 4.9E24/M_sun;
	double M_earth = 6.0E24/M_sun;
	double M_mars = 6.6E23/M_sun;  
	double M_jupiter = 1.9E27/M_sun;
	double M_saturn = 5.5E26/M_sun; 
	double M_uranus = 8.8E25/M_sun;
	double M_neptun = 1.03E26/M_sun; 
	double M_pluto = 1.31E22/M_sun;     

	// Intialize planets 	

	planet sun(1.0, 0, 0, 0, 0, "Sun"); 
	planet mercury(M_mercury,-3.921538113765745E-01,-5.487782323669228E-02,-1.585262467570905E-03*year,-2.659352789459454E-02*year,"Mercury");
	planet venus(M_venus,-5.262376668671904E-01, 4.909161410927460E-01, -1.375950023480965E-02*year, -1.499877247340715E-02*year, "Venus"); 
	planet earth(M_earth, 9.734E-01, 2.4198E-01, -4.349E-03*year, 1.665E-02*year, "Earth"); 
	planet mars(M_mars, -1.5158E+00, 6.9015E-01, -5.2319E-03*year, -1.1556E-02*year, "Mars"); 
	planet jupiter(M_jupiter, -4.623E+00, -2.860E+00, 3.881E-03*year, -6.058E-03*year, "Jupiter"); 
	planet saturn(M_saturn,-4.054191284158095E-01,-1.004694419069413E+01,5.268688859602635E-03*year,-2.426734651712632E-04*year,"Saturn"); 
	planet uranus(M_uranus,1.787729807481248E+01,8.776369655211150E+00,-1.762081973967272E-03*year,3.347240820361182E-03*year,"Uranus"); 
	planet neptun(M_neptun,2.860478537542545E+01,-8.851521431448582E+00,9.066379590607452E-04*year,3.017237536736915E-03*year,"Neptun"); 
	planet pluto(M_pluto,1.051657650590531E+01,-3.171589804117349E+01,3.055852032837181E-03*year,3.423298915677500E-04*year,"Pluto"); 
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

	solver Binary; 
	Binary.addPlanet(sun); 
	Binary.addPlanet(earth);

	clock_t FE_start, FE_finish; 
	FE_start = clock();  
	Binary.ForwardEuler(n, t);  
	FE_finish = clock(); 

	clock_t VV_start, VV_finish; 
	VV_start = clock(); 
	Binary.VelocityVerlet(n, t); 
	VV_finish = clock(); 

	// Three-body system (earth, sun and jupiter) 

	solver ThreeBody; 
	ThreeBody.addPlanet(sun);
	ThreeBody.addPlanet(earth); 
	ThreeBody.addPlanet(jupiter);  
	ThreeBody.VelocityVerlet(n, t); 

	// Solar system (all planets) 

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


	cout << "----------------------" << endl; 

	double FE_time = (FE_finish - FE_start)/((double) CLOCKS_PER_SEC); 
	double VV_time = (VV_finish - VV_start)/((double) CLOCKS_PER_SEC); 	

	cout << "Time with FE: " << FE_time << endl; 
	cout << "Time with VV: " << VV_time << endl; 

}
