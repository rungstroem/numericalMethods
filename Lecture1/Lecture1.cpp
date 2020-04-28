#include <iostream>
#include "nr3.h"
#include "ludcmp.h"
#include "utilities.h"

using namespace std;
	
int main() {



	// Exercise 1:
	// Solve A x = b using LU decomposition, and print the result.

	MatDoub A(3,3);
	A[0][0] = 1.0;	A[0][1] = 2.0;	A[0][2] = 3.0;
	A[1][0] = 2.0;	A[1][1] = -4.0;	A[1][2] = 6.0;
	A[2][0] = 3.0;	A[2][1] = -9.0;	A[2][2] = -3.0;

	VecDoub b(3);
	b[0] = 5.0;
	b[1] = 18.0;
	b[2] = 6.0;


	// evaluate x

	// Create LU-decomt object
	LUdcmp lu(A);

	//Output vector x
	VecDoub_O x(3);

	//Call LU-decomp solve function
	lu.solve(b,x);
	
	// print x
	util::print(lu.lu, "LU Matrix");
	util::print(A,"Matrix A");
	util::print(x, "Vector x");
	util::print(b, "Vector b");

	//Decompose LU matrix in Alpha & Beta matrices
	MatDoub LU = lu.lu;
    MatDoub a(3,3);     //Placeholder for resulterende A-matrix fra L*U
	
    //Beta Matricen er den øvre matrix - ie. U del ad LU-decmop.
    MatDoub Beta(3,3);

	Beta[0][0] = LU[0][0];	Beta[0][1] = LU[0][1];	Beta[0][2] = LU[0][2];
	Beta[1][0] = 0;		    Beta[1][1] = LU[1][1];	Beta[1][2] = LU[1][2];
	Beta[2][0] = 0;		    Beta[2][1] = 0;		    Beta[2][2] = LU[2][2];
    
    //Alpha matricen er den nedre matrix - ie. L del af LU-decomp.
    //Alpha matricen skal have 1'er på diagonalen
    MatDoub Alpha(3,3);
    
	Alpha[0][0] = 1;	    Alpha[0][1] = 0;	    Alpha[0][2] = 0;
	Alpha[1][0] = LU[1][0];	Alpha[1][1] = 1;	    Alpha[1][2] = 0;
	Alpha[2][0] = LU[2][0];	Alpha[2][1] = LU[2][1];	Alpha[2][2] = 1;

	//util::print(Beta, "Beta Matrix");
	//util::print(Alpha, "Alpha Matrix");
	

	//Matrix[ROW][COL]	->	Matrix[i][j]

	//MatDoub Alpha(3,3); MatDoub Beta(3,3);

    /*
	Alpha[0][0] = 1;	Alpha[0][1] = 2;	Alpha[0][2] = 3;	
	Alpha[1][0] = 4;	Alpha[1][1] = 5;	Alpha[1][2] = 6;
	Alpha[2][0] = 7;	Alpha[2][1] = 8;	Alpha[2][2] = 9;

	Beta[0][0] = 2;		Beta[0][1] = 2;		Beta[0][2] = 2;
	Beta[1][0] = 3;		Beta[1][1] = 3;		Beta[1][2] = 3;
	Beta[2][0] = 4;		Beta[2][1] = 4;		Beta[2][2] = 4;
    
	//a[0][0] = Alpha[0][0]*Beta[0][0] + Alpha[0][1]*Beta[1][0];
	//a[0][1] = Alpha[0][0]*Beta[0][1] + Alpha[0][1]*Beta[1][1];

	//a[1][0] = Alpha[1][0]*Beta[0][0] + Alpha[1][1]*Beta[1][0];
	//a[1][1] = Alpha[1][0]*Beta[0][1] + Alpha[1][1]*Beta[1][1];

	//a[i0][j0] = Alpha[i0][j0]*Beta[i0][j0] + Alpha[i0][j1]*Beta[i1][j0];
	//a[i0][j1] = Alpha[i0][j0]*Beta[i0][j1] + Alpha[i0][j1]*Beta[i1][j1];
	
	//a[1][0] = Alpha[1][0]*Beta[0][0] + Alpha[1][1]*Beta[1][0];
	//a[1][1] = Alpha[1][0]*Beta[0][1] + Alpha[1][1]*Beta[1][1];
	

	//Working matrix multiplication - yet no dynamic size!!
	int row = 3;	int col = 3;	

	for(int i = 0;i<row;i++){

		for(int k = col;k<3;k++){
		for(int j = col;j<3;j++){
			a[i][k] = a[i][k] + Alpha[i][j]*Beta[j][i];
		}
		}
		
	}
    */

    // Remember!! A = L*U - ei. A = Alpha*Beta
    a =operator*(Alpha, Beta);
	//Print the resulting matrix from multiplication above
	util::print(a, "Matrix A from L*U");	
	return 0;
}
