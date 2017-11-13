#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <string>
#include <iomanip>
#include <armadillo>
#include "time.h"
#include "lib.h"

using namespace std;
using namespace arma; 

inline int periodic(int i, int limit, int add){
	return (i+limit+add)%limit; 
}

void initialize(int, double, mat, double&, double&);
void Metropolis(int, long&, mat&, double&, double&, vec&); 

int main(int argc, char *argv[]){ 

	long idum; 
	int n_spins = atoi(argv[1]); 
	int mcs = atoi(argv[2]);
	mat spin_matrix = ones<mat>(n_spins, n_spins);  
	double initial_temp, final_temp, temp_step, E, M;
	vec w = zeros<vec>(17); 
	vec average = zeros<vec>(5); 

	initial_temp = final_temp = temp_step = 1.0; 
	idum = -1; 
	for(double temp = initial_temp; temp <= final_temp; temp+=temp_step){
		E = M = 0; 	
		for( int de = -8; de <= 8; de++) w(de+8) = 0; 
		for( int de = -8; de <= 8; de+=4) w(de+8) = exp(-de/temp); 

		initialize(n_spins, temp, spin_matrix, E, M);		

		for( int cycles = 1; cycles <= mcs; cycles++){
			Metropolis(n_spins, idum, spin_matrix, E, M, w); 
			average(0) += E; average(1) += E*E; average(2) += M; average(3) += M*M; average(4) += fabs(M); 
		}

	}	

	for( int i = 0; i<5; i++) average(i) = average(i)/((double) mcs); 

	double C_v = average(1) - average(0)*average(0); 
	double chi = average(2) - average(3)*average(3); 

	cout << "--------------" << endl; 
	cout << "Output: " << endl;   
	cout << "Energy: " << average(0) << endl; 
	cout << "Magnetization: " << average(4) << endl; 
	cout << "Specific heat: " << C_v << endl; 
	cout << "Suceptibility: " << chi << endl; 


	return 0; 

}


void initialize(int n_spins, double temperature, mat spin_matrix, double &E, double &M){

	for(int i = 0; i<n_spins; i++){
		for(int j = 0; j<n_spins; j++){
			M += (double) spin_matrix(i,j); 
			E -= (double) spin_matrix(i,j)*( spin_matrix(periodic(j,n_spins,-1), i) + spin_matrix(j, periodic(i,n_spins,-1) )  ); 
		}
	}

}


void Metropolis(int n_spins, long& idum, mat& spin_matrix, double& E, double& M, vec& w){

	for(int y = 0; y<n_spins; y++){
		for(int x = 0; x<n_spins; x++){
			int ix = (int) (ran1(&idum)*(double)n_spins); 
			int iy = (int) (ran1(&idum)*(double)n_spins); 

			int deltaE = 2*spin_matrix(iy,ix)*(  spin_matrix(iy, periodic(ix,n_spins,-1)) + spin_matrix(periodic(iy,n_spins,-1),ix)	
										 										+	 spin_matrix(iy, periodic(ix,n_spins, 1)) + spin_matrix(periodic(iy,n_spins, 1),ix) );

			if( ran1(&idum) <= w(deltaE + 8) ){
				spin_matrix(iy,ix) *= -1; 
				M += (double) 2*spin_matrix(iy, ix); 
				E += (double) deltaE; 
			}

		}
	}			

}



