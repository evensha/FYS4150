#include "planet.h"

planet::planet(){  
	// Initialize planet with no specifications 

	mass = 0.0; 
	position[0] = 0.0; 
	position[1] = 0.0; 
	position[2] = 0.0; 
	velocity[0] = 0.0; 
	velocity[1] = 0.0; 
	velocity[2] = 0.0; 

	name = "Planet"; 	

}

planet::planet(double M, double x, double y, double vx, double vy){ 
	// Initialize planet with mass and two dimensions   

	mass = M; 
	position[0] = x; 
	position[1] = y;  
	position[2] = 0.0; 
	velocity[0] = vx; 
	velocity[1] = vy; 
	velocity[2] = 0.0; 

	name = "Planet"; 

}


planet::planet(double M, double x, double y, double vx, double vy, string planet){  
	// Initialize planet with mass, two dimensions and name of the planet    

	mass = M; 
	position[0] = x; 
	position[1] = y;  
	position[2] = 0.0; 
	velocity[0] = vx; 
	velocity[1] = vy; 
	velocity[2] = 0.0; 

	name = planet; 

}


planet::planet(double M, double x, double y, double z, double vx, double vy, double vz, string planet){ 
	// Initialize planet with mass, three dimensions and name of the planet    

	mass = M; 
	position[0] = x; 
	position[1] = y;  
	position[2] = z; 
	velocity[0] = vx; 
	velocity[1] = vy; 
	velocity[2] = vz; 

	name = planet; 

}


double planet::Distance(planet otherPlanet){
	// Calculate the distance between this planet and another planet 

	double x1,y1,z1,x2,y2,z2,Dx,Dy,Dz; 
	
	x1 = this->position[0]; 
	y1 = this->position[1]; 		
	z1 = this->position[2]; 

	x2 = otherPlanet.position[0]; 
	y2 = otherPlanet.position[1]; 
	z2 = otherPlanet.position[2]; 	

	Dx = x1 - x2; 
	Dy = y1 - y2; 
	Dz = z1 - z2; 

	return sqrt(Dx*Dx + Dy*Dy + Dz*Dz); 

}


double planet::KineticEnergy(){ 
	// Calculate the planets kinetic energy 

	double v_x = this->velocity[0]; double v_y = this->velocity[1]; double v_z = this->velocity[2]; 
	double m = this->mass; 		
	double E_k = 0.5*m*(v_x*v_x + v_y*v_y + v_z*v_z); 

	return E_k; 

}


double planet::PotentialEnergy(planet otherPlanet){
	// Calculate the planets kinetic energy 

	double m = this->mass; 
	double M = otherPlanet.mass; 
	double r = this->Distance(otherPlanet);  

	double E_p = - 4*M_PI*M_PI*M*m/r; 
	 	
	return E_p; 


}

double planet::xMomentum(){
	// Calculate the planets momentum in x-direction 

	double v_x = this->velocity[0]; 
	double m = this->mass; 

	return v_x*m; 
}

double planet::yMomentum(){
	// Calculate the planets momentum in y-direction 

	double v_y = this->velocity[1]; 
	double m = this->mass; 

	return v_y*m; 
}

double planet::AngularMomentum(){
	// Calculate the planets angular momentum 

	double x = this->position[0]; double y = this->position[1]; double z = this->position[2]; 
	double v_x = this->velocity[0]; double v_y = this->velocity[1]; double v_z = this->velocity[2]; 

	double L_x = y*v_z - z*v_y; 
	double L_y = -(x*v_z - z*v_x);  
	double L_z = x*v_y - y*v_x; 
	double L_tot = sqrt(L_x*L_x + L_y*L_y + L_z*L_z); 

	return L_tot;  

}




