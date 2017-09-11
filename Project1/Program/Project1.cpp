#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <iomanip>

using namespace std;

ofstream ofile; 

int main(int argc, char* argv[]){

	char *outfilename;
	outfilename = argv[1];
	ofile.open(outfilename); 

	int n = atof(argv[2]); 

	double h = 1.0/((double n)+1); 

	a = double matr[n]; 	

	for(int i=0; i<n; i++){ 
		a[i] = 2;  
	} 
	cout << a << endl; 
	ofile.close(); 
	return 0; 
}  
