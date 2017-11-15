#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <string>
#include <iomanip>
#include <armadillo>
#include "time.h"
#include "lib.h"
//#include "mpi.h"

using namespace std;
using namespace arma; 

ofstream ofile; 

inline int periodic(int i, int limit, int add){
	return (i+limit+add) % (limit); 
}


void initialize(int, double, mat, double&, double&);
void Metropolis(int, long&, mat&, double&, double&, vec); 
void output(int, int, double, vec); 

int main(int argc, char *argv[]){ 

	// Defining quantities and the spin matrix 

	long idum; 
	int n_spins = atoi(argv[1]); 
	int mcs = atoi(argv[2]); 
	int rand_spins = atoi(argv[3]); 
	mat spin_matrix; 
	string spins = "ordered";   
	idum = -1; 

	double initial_temp, final_temp, temp_step, E, M;
	vec w = zeros<vec>(33); 
	vec average = zeros<vec>(5); 
	vec total_average = zeros<vec>(5); 
	int my_rank = 0; 
	/*
	// MPI initialization 
	int my_rank, numprocs; 
	MPI_Init(&argc, &argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs); 
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
	if(my_rank == 0 && argc <= 1){
		cout << "Bad usage: " << argv[0] << " read output file" << endl; exit(1); 
	}
	


	int no_intervalls = mcs/numprocs; 
	int myloop_begin = my_rank*no_intervalls +1; 
	int myloop_end = (my_rank+1)*no_intervalls; 
	if( (my_rank == numprocs-1) && (myloop_end < mcs ) ) myloop_end = mcs; 

	// Broadcast to all nodes common variables 

	MPI_Bcast(&n_spins, 1, MPI_INT, 0, MPI_COMM_WORLD); 
	MPI_Bcast(&initial_temp, 1, MPI_INT, 0, MPI_COMM_WORLD); 
	MPI_Bcast(&final_temp, 1, MPI_INT, 0, MPI_COMM_WORLD); 
	MPI_Bcast(&temp_step, 1, MPI_INT, 0, MPI_COMM_WORLD); 

	idum = -1-my_rank; 
	double TimeStart, TimeEnd, TotalTime; 
	TimeStart = MPI_Wtime(); 
	*/
	initial_temp = 1.0; final_temp = 2.4; temp_step = 1.4; 
	for(double temp = initial_temp; temp <= final_temp; temp+=temp_step){

		E = M = 0; 	
		spin_matrix = ones<mat>(n_spins, n_spins);
		if( rand_spins  == 1 ){
			spins = "random"; 
			for(int x = 0; x<n_spins; x++){
				for(int y = 0; y<n_spins; y++){
					double s = ran1(&idum); 
					if( s < 0.5 ) spin_matrix(x,y) = -1; 
					if( s > 0.5 ) spin_matrix(x,y) = 1; 
				}
			}
		}

		for( int de = -16; de <= 16; de++) w(de+16) = 0; 
		for( int de = -16; de <= 16; de+=8) w(de+16) = exp(-de/temp); 
		for( int i = 0; i<5; i++) average(i) = 0; 

		initialize(n_spins, temp, spin_matrix, E, M);	
		/*
		cout << "--------------------" << endl; 	
		cout << "Initial values:" << endl; 
		cout << "Energy: " << E << endl; 
		cout << "Magnetization: " << M << endl; 
		*/
		ostringstream os; 
		os << "Output/Output_" << spins << "_" << temp << ".txt"; 
		string outfile = os.str(); 
		ofile.open(outfile.c_str()); 
		ofile << "Time" << setw(20) << "Energy" << setw(20) << "Magnetization" << endl;  

		for( int cycles = 1; cycles <= mcs; cycles++){
			Metropolis(n_spins, idum, spin_matrix, E, M, w); 
			average(0) += E; average(1) += E*E; average(2) += M; average(3) += M*M; average(4) += fabs(M); 
			//output(n_spins, mcs, temp, total_average); 
			ofile << cycles << setw(20) << average(0)/((double) cycles) << setw(20) << average(4)/((double) cycles) << endl; 
		}
		/*	
		// Find total average
		for( int i = 0; i<5; i++){
			MPI_Reduce(&average(i), &total_average(i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
		}
		
		if( my_rank == 0 ){
			output(n_spins, mcs, temp, total_average); 
		}
		*/
		ofile.close(); 
		//TimeEnd = MPI_Wtime(); 
		//TotalTime = TimeEnd - TimeStart; 
		
		if( my_rank == 0 ){

			cout << "--------------------" << endl; 
			//cout << "Time = " << TotalTime << " on number of processors: " << numprocs << endl; 
			cout << "Temperature = " << temp << endl; 	

			for( int i = 0; i<5; i++) average(i) = average(i)/((double) mcs); 	
			double C_v = average(1) - average(0)*average(0); 
			double chi = average(3) - average(2)*average(2); 

			cout << "Final values (averages): " << endl;   
			cout << "Energy: " << average(0) << endl; 
			//cout << "Magnetization (mean): " << average(3) << endl; 
			cout << "Magnetization (abs): " << average(4) << endl; 
			cout << "Specific heat: " << C_v << endl; 
			cout << "Suceptibility: " << chi << endl; 
			cout << "--------------------" << endl; 

		}
		
	}	

	//MPI_Finalize(); 
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


void Metropolis(int n_spins, long& idum, mat& spin_matrix, double& E, double& M, vec w){

	for(int y = 0; y<n_spins; y++){
		for(int x = 0; x<n_spins; x++){
			int ix = (int) (ran1(&idum)*(double)n_spins); 
			int iy = (int) (ran1(&idum)*(double)n_spins); 
			int deltaE = 2*spin_matrix(iy,ix)*(  spin_matrix(iy, periodic(ix,n_spins,-1)) + spin_matrix(periodic(iy,n_spins,-1),ix)	
										 										+	 spin_matrix(iy, periodic(ix,n_spins, 1)) + spin_matrix(periodic(iy,n_spins, 1),ix) );
			if( ran1(&idum) <= w(deltaE + 16) ){
				spin_matrix(iy,ix) *= -1; 
				M += (double) 2*spin_matrix(iy, ix); 
				E += (double) deltaE; 
			}

		}
	}			

}


void output(int n_spins, int mcs, double temperature, vec total_average){

	double norm = 1/((double) (mcs)); 
	double Etotal_average = total_average(0)*norm; 
	double E2total_average = total_average(1)*norm; 
	double Mtotal_average = total_average(2)*norm; 
	double M2total_average = total_average(3)*norm; 
	double Mabstotal_average = total_average(4)*norm; 

	double Evariance = (E2total_average - Etotal_average*Etotal_average)/n_spins/n_spins; 
	double Mvariance = (M2total_average - Mtotal_average*Mtotal_average)/n_spins/n_spins; 

	ofile << setiosflags(ios::showpoint | ios::uppercase); 
	ofile << setw(15) << setprecision(8) << temperature; 
	ofile << setw(15) << setprecision(8) << Etotal_average/n_spins/n_spins;
	ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature; 
	ofile << setw(15) << setprecision(8) << Mtotal_average/n_spins/n_spins; 
	ofile << setw(15) << setprecision(8) << Mvariance/temperature; 
	ofile << setw(15) << setprecision(8) << Mabstotal_average/n_spins/n_spins << endl; 

}







