#ifndef SOLVER_H
#define SOLVER_H 
#include "planet.h"
#include <fstream>
#include <vector> 


using std::vector; 

class solver{ 

public: 
	friend class planet; 

	// Properties 
	double mass, G;
	double radius; 
	int total_planets; 
	vector<planet> all_planets; 

	// Initializers 
	solver(); 
	solver(double r); 

	// Functions 
	void addPlanet(planet newplanet); 
	void ForwardEuler(int integration_points, double time); 
	void VelocityVerlet(int integration_points, double time); 
	void GravitationalForce(planet &Planet, planet &other, double &F_x, double &F_y); 
	void PrintPositions(); 
	void PrintNames(); 

}; 

#endif // SOLVER_H
