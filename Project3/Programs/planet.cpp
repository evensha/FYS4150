/*
#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <string>
#include <iomanip>
#include <armadillo>
#include "time.h"
*/
#include "planet.h"

//using namespace std;
//using namespace arma; 

planet::planet(){

	mass = 0.0; 
	position[0] = 0.0; 
	position[1] = 0.0; 
	velocity[0] = 0.0; 
	velocity[1] = 0.0; 
	kinetic = 0.0; 
	potential = 0.0; 

	name = "some_planet"; 	

}

planet::planet(double M, double x, double y, double vx, double vy){ 

	mass = M; 
	position[0] = x; 
	position[1] = y;  
	velocity[0] = vx; 
	velocity[1] = vy; 
	kinetic = 0; 
	potential = 0; 

	name = "some_planet"; 

}


planet::planet(double M, double x, double y, double vx, double vy, string planet){ 

	mass = M; 
	position[0] = x; 
	position[1] = y;  
	velocity[0] = vx; 
	velocity[1] = vy; 
	kinetic = 0; 
	potential = 0; 

	name = planet; 

}


double planet::Distance(planet otherPlanet){

	double x1,y1,x2,y2,Dx,Dy; 
	
	x1 = this->position[0]; 
	y1 = this->position[1]; 		

	x2 = otherPlanet.position[0]; 
	y2 = otherPlanet.position[1]; 
	
	Dx = x1 - x2; 
	Dy = y1 - y2; 

	return sqrt(Dx*Dx + Dy*Dy); 

}


double planet::GravitationalForce(planet otherPlanet, double Gconst){
	
	double r = this->Distance(otherPlanet); 

	if(r!=0) return Gconst*this->mass*otherPlanet.mass/(r*r); 

	else return 0; 

}


double planet::Acceleration(planet otherPlanet, double Gconst){

	double r = this->Distance(otherPlanet); 

	if(r!=0) return this->GravitationalForce(otherPlanet, Gconst)/(this->mass*r); 

	else return 0; 

}


double planet::KineticEnergy(){ 

	double v_x = this->velocity[0]; double v_y = this->velocity[1]; 
	double m = this->mass; 		
	double E_k = 0.5*m*(v_x*v_x + v_y*v_y); 

	return E_k; 

}


double planet::PotentialEnergy(double Gconst){

	double m = this->mass; double x = this->position[0]; double y = this->position[1]; 
	double r = sqrt(x*x + y*y); 

	double E_p = - Gconst*m/r; 
	 	
	return E_p; 

}


double planet::AngularMomentum(){

	double x = this->position[0]; double y = this->position[1]; 
	double v_x = this->position[0]; double v_y = this->position[1]; 

	double L_z = x*v_y - y*v_x; 
	double L_tot = sqrt(L_z*L_z); 

	return L_tot;  

}




