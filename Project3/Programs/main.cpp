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

	planet earth(0.000003, 1.0, 0.0, 0.0, v_start*2*M_PI); 
	//planet sun(1.0, 0, 0, 0, 0); 

	solver EarthSun_FE(1.0);
	EarthSun_FE.addPlanet(earth); 
	EarthSun_FE.ForwardEuler(n, t); 

	solver EarthSun_VV(1.0); 
	EarthSun_VV.addPlanet(earth); 
	EarthSun_VV.VelocityVerlet(n, t); 

}
