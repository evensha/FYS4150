#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma; 

double offdiag( mat A, int *p, int *q, int n ); // function that check value of the off-diagonal elements 
void Jacobi_rotation( mat& A, mat& R, int k, int l, int n ); // function that does the rotation 

void do_Jacobi(mat& A, mat& R, vec& lambda, int n); // function that runs the algorithm

