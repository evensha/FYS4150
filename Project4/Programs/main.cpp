#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <string>
#include <iomanip>
#include <armadillo>
#include "time.h"
#include "lib.h"
#include "mpi.h"

using namespace std;
using namespace arma; 

ofstream ofile, ofile1; 

inline int periodic(int i, int limit, int add){
	return (i+limit+add) % (limit); 
}


void initialize(int, double, mat, double&, double&);
void Metropolis(int, long&, mat&, double&, double&, int&, vec); 
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
	//int my_rank = 0; 
	
	// MPI initialization 
	int my_rank, numprocs; 
	MPI_Init(&argc, &argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs); 
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
	if(my_rank == 0 && argc <= 1){
		cout << "Bad usage: " << argv[0] << " read output file" << endl; exit(1); 
	}
	if( my_rank == 0 && argc > 1 ){
		ostringstream os1; 
		os1 << "Output/Output_L_" << n_spins << ".txt"; 
		string outfile1 = os1.str(); 
		ofile1.open(outfile1.c_str()); 
	}

	if( n_spins == 20){
		 initial_temp = 1.0; final_temp = 2.4; temp_step = 1.4; 
	}
	else{ 
		initial_temp = 2.0; final_temp = 2.3; temp_step = 0.05; 
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

	int acc_confs; 
	for(double temp = initial_temp; temp <= final_temp; temp+=temp_step){

		E = M = acc_confs = 0; 	
		spin_matrix = ones<mat>(n_spins, n_spins);
		if( rand_spins  == 1 ){
			spins = "random"; 
			for(int x = 0; x<n_spins; x++){
				for(int y = 0; y<n_spins; y++){
					double s = ran1(&idum); 
					if( s < 0.5 ) spin_matrix(x,y) = -1; 
					else spin_matrix(x,y) = 1; 
					//cout << spin_matrix(x,y) << endl; 
				}
			}
		}

		for( int de = -16; de <= 16; de++) w(de+16) = 0; 
		for( int de = -16; de <= 16; de+=8) w(de+16) = exp(-de/temp); 
		for( int i = 0; i<5; i++) average(i) = 0; 
		for( int i = 0; i<5; i++) total_average(i) = 0; 

		initialize(n_spins, temp, spin_matrix, E, M);	
		/*
		cout << "--------------------" << endl; 	
		cout << "Initial values:" << endl; 
		cout << "Energy: " << E << endl; 
		cout << "Magnetization: " << M << endl; 
		*/
		if( n_spins == 20){
			ostringstream os; 
			os << "Output/Output_" << spins << "_" << temp << ".txt"; 
			string outfile = os.str(); 
			ofile.open(outfile.c_str()); 
			ofile << "MCs" << setw(20) << "E" << setw(20) << "meanE" << setw(20) << "M" << setw(20) << "meanM" << setw(20) << 
			"Accepted configs" << endl;  
		}
	
		for( int cycles = myloop_begin; cycles <= myloop_end; cycles++){
			Metropolis(n_spins, idum, spin_matrix, E, M, acc_confs, w); 
			average(0) += E; average(1) += E*E; average(2) += M; average(3) += M*M; average(4) += fabs(M); 
			if( n_spins == 20){
			ofile << cycles << setw(20) << E << setw(20) << average(0)/((double) cycles) << setw(20) << 
																		 M << setw(20) << average(4)/((double) cycles) << setw(20) << acc_confs << endl; 
			}
		}

		// Find total average
		for( int i = 0; i<5; i++){
			MPI_Reduce(&average(i), &total_average(i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
		}
		
		if( my_rank == 0 ){
			output(n_spins, mcs, temp, total_average); 
		}
		
		//output(n_spins, mcs, temp, average); 
		if(n_spins == 20) ofile.close(); 
		TimeEnd = MPI_Wtime(); 
		TotalTime = TimeEnd - TimeStart; 
		
		if( my_rank == 0 ){

			cout << "--------------------" << endl; 
			cout << "Time = " << TotalTime << " on number of processors: " << numprocs << endl; 
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

	MPI_Finalize(); 
	ofile1.close(); 
	return 0; 

}


void initialize(int n_spins, double temperature, mat spin_matrix, double &E, double &M){

	for(int i = 0; i<n_spins; i++){
		for(int j = 0; j<n_spins; j++){
			M += (double) spin_matrix(i,j); 
			E -= (double) spin_matrix(i,j)*( spin_matrix(periodic(i,n_spins,-1), j) + spin_matrix(i, periodic(j,n_spins,-1) )  ); 
		}
	}

}


void Metropolis(int n_spins, long& idum, mat& spin_matrix, double& E, double& M, int& a, vec w){

	for(int y = 0; y<n_spins; y++){
		for(int x = 0; x<n_spins; x++){
			int ix = (int) (ran1(&idum)*(double)n_spins); 
			int iy = (int) (ran1(&idum)*(double)n_spins); 
			int deltaE = 2*spin_matrix(iy,ix)*(  spin_matrix(iy, periodic(ix,n_spins,-1)) + spin_matrix(periodic(iy,n_spins,-1),ix)	
										 										+	 spin_matrix(iy, periodic(ix,n_spins, 1)) + spin_matrix(periodic(iy,n_spins, 1),ix) );
			//if( deltaE == 32 ) cout << "What!!!" << endl; 
			if( ran1(&idum) <= w(deltaE + 16) ){
				a++; 
				spin_matrix(iy,ix) *= -1; 
				M += (double) 2*spin_matrix(iy, ix); 
				E += (double) deltaE; 
			}

		}
	}			

}


void output(int n_spins, int mcs, double temperature, vec total_average){
	//cout << "Inside output function!!" << endl; 
	double norm = 1/((double) (mcs)); 
	double Etotal_average = total_average(0)*norm; 
	double E2total_average = total_average(1)*norm; 
	double Mtotal_average = total_average(2)*norm; 
	double M2total_average = total_average(3)*norm; 
	double Mabstotal_average = total_average(4)*norm; 

	double Evariance = (E2total_average - Etotal_average*Etotal_average); ///n_spins/n_spins; 
	double Mvariance = (M2total_average - Mtotal_average*Mtotal_average); ///n_spins/n_spins; 

	ofile1 << setiosflags(ios::showpoint | ios::uppercase); 
	ofile1 << setw(15) << setprecision(8) << temperature; 
	ofile1 << setw(15) << setprecision(8) << Etotal_average; ///n_spins/n_spins;
	ofile1 << setw(15) << setprecision(8) << Evariance/temperature/temperature; 
	ofile1 << setw(15) << setprecision(8) << Mtotal_average; ///n_spins/n_spins; 
	ofile1 << setw(15) << setprecision(8) << Mvariance/temperature; 
	ofile1 << setw(15) << setprecision(8) << Mabstotal_average << endl; ///n_spins/n_spins << endl; 

}







