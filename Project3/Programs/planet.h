#ifndef PLANET_H 
#define PLANET_H
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <string>
using std::vector; 
using std::string; 
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
	string name; 

	// Initializers 
	planet(); 
	planet(double M, double x, double y, double vx, double vy); 
	planet(double M, double x, double y, double vx, double vy, string planet); 

	// Functions 	
	double Distance(planet otherPlanet); 
	double GravitationalForce(planet otherPlanet, double Gconst);
	double Acceleration(planet otherPlanet, double Gconst); 
	double KineticEnergy(); 
	double PotentialEnergy(double Gconst); 
	double AngularMomentum(); 


};

#endif  // PLANET_H
