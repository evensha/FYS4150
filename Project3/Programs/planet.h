#ifndef PLANET_H 
#define PLANET_H
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
using std::vector; 
/*
#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <string>
#include <iomanip>
*/
//#include <armadillo>
//#include "time.h"

class planet{ 

public: 

	// Properties 
	double mass; 
	double position[2]; 
	double velocity[2]; 
	double kinetic; 
	double potential; 

	// Initializers 
	planet(); 
	planet(double M, double x, double y, double vx, double vy); 

	// Functions 	
	double distance(planet otherPlanet); 
	double GravitationalForce(planet otherPlanet, double Gconst);
	double acceleration(planet otherPlanet, double Gconst); 
	double KineticEnergy(); 
	double PotentialEnergy(planet &otherPlanet, double Gconst); 

};

#endif  // PLANET_H
