#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace std;

ofstream ofile; 

int main(int argc, char* argv[]){

	char *outfilename;
	outfilename = argv[1];
	ofile.open(outfilename); 

	int n = atof(argv[2]); 

	//double h = 1.0/((double n)+1); 
	/*
	double A[n][n]; 
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(i==j){A[i][j] = 2;} 
			else{A[i][j] = 0;} 
		}
	}

	cout << A[1][1] << "," << A[n-1][n-1] << "," << A[3][4] << endl;  
	*/
	ofile.close(); 
	return 0; 
}  
