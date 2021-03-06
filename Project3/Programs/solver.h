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
	double mass, G, beta;
	int RelCorr; 
	int total_planets; 
	vector<planet> all_planets; 

	// Initializers 
	solver(); 
	solver(double b);
	solver(int withGR);  

	// Functions 
	void addPlanet(planet newplanet); 
	void ForwardEuler(int integration_points, double time, int withOutput); 
	void VelocityVerlet(int integration_points, double final_time, int withOutput); 
	void GravitationalForce(planet &Planet, planet &other, double &F_x, double &F_y, double &F_z, double beta, int RelCorr); 
	void PrintPositions(); 
	void PrintNames(); 

}; 

#endif // SOLVER_H
