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

	int n = atoi(argv[1]);
	double t = atof(argv[2]); 

	double v_start = atof(argv[3]);   

	//planet earth(0.000003, 9.734E-01, 2.4198E-01, -4.349E-03*365, 1.665E-02*365); 
	planet earth(0.000003, 1.0, 0.0, 0.0, v_start*2*M_PI); 	
	planet sun(1.0, 0, 0, 0, 0); 

	cout << "----------------------" << endl; 
	cout << "Initial values: " << endl; 
	cout << "Kinetic energy: " << earth.KineticEnergy() << endl; 
	cout << "Potential energy: " << earth.PotentialEnergy(4*M_PI*M_PI) << endl; 
	cout << "Angular momentum: " << earth.AngularMomentum() << endl; 

	clock_t FE_start, FE_finish; 
	FE_start = clock(); 

	solver EarthSun_FE(1.0); 
	EarthSun_FE.addPlanet(earth); 
	EarthSun_FE.ForwardEuler(n, t);  

	FE_finish = clock(); 

	clock_t VV_start, VV_finish; 
	VV_start = clock(); 

	solver EarthSun_VV(1.0); 
	EarthSun_VV.addPlanet(earth); 
	EarthSun_VV.addPlanet(sun); 
	EarthSun_VV.VelocityVerlet(n, t); 

	VV_finish = clock(); 

	cout << "----------------------" << endl; 

	double FE_time = (FE_finish - FE_start)/((double) CLOCKS_PER_SEC); 
	double VV_time = (VV_finish - VV_start)/((double) CLOCKS_PER_SEC); 	

	cout << "Time with FE: " << FE_time << endl; 
	cout << "Time with VV: " << VV_time << endl; 

}
