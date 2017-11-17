#ifndef PLANET_H 
#define PLANET_H
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <string>
using std::vector; 
using std::string; 

class planet{ 

public: 

	// Properties 
	double mass; 
	double position[3]; 
	double velocity[3]; 
	string name; 

	// Initializers 
	planet(); 
	planet(double M, double x, double y, double vx, double vy); 
	planet(double M, double x, double y, double vx, double vy, string planet); 
	planet(double M, double x, double y, double z, double vx, double vy, double vz, string planet); 
	
	// Functions 	
	double Distance(planet otherPlanet); 
	double KineticEnergy(); 
	double PotentialEnergy(planet otherPlanet); 
	double xMomentum();
	double yMomentum();  
	double AngularMomentum(); 


};

#endif  // PLANET_H
