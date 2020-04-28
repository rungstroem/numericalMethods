#include <iostream>
#include <fstream>
#include "nr3.h"
#include "ludcmp.h"
#include "cholesky.h"
#include "utilities.h"
#include "svd.h"
#include <math.h>
using namespace std;

void pont_data(){

VecDoub xPont(40); VecDoub yPont(40);
ifstream Pont("PontiusData.dat");
for(int i = 0; i < 40; i++) {
	Pont >> yPont[i];
	Pont >> xPont[i];
}

// Normal Equation for solving problem
// (AT*A)*a = AT*b
// A design matrix,
// AT transpose of A
// a the vector to solve fore - typically named x
// b vector with result from equation A*x=b

// Skal bruges til at få std senere til at passe - Divider matrix A elementer
// og vektor b elementer
Doub std_dev_pont = 0.0002051;

// Creating matrix A from data, object xPont
MatDoub A(40,3);
for(int i=0;i<40;i++){
	A[i][0] = 1/std_dev_pont;
	A[i][1] = xPont[i]/std_dev_pont;
	A[i][2] = xPont[i]*xPont[i]/std_dev_pont;
}

// Creating vector b from data, object yPont
VecDoub b(40);
for(int i=0;i<40;i++){
	b[i] = yPont[i]/std_dev_pont;
}

// Transposing matrix A
MatDoub AT= util::Transpose(A);

// From normal equation
// AT*A = C
// Capital C is the resulting matrix
// In book capital C is alpha
//
// AT*b = c
// lowercase c is the resulting vector
// In book lowercase c is beta

// Creating matrix C
MatDoub C = operator*(AT, A);
// Creating vector c
VecDoub c = operator*(AT, b);

// LU decomposition solve of linear equation problem
LUdcmp lu(C);
VecDoub_O x(3);		//a in (AT*A)*a = AT*b	C*a=c

lu.solve(c, x);

// Cholesky solve of the linear equation problem
Cholesky C_cho(C);
VecDoub_O x_cho(3);	//a in (AT*A)*a = AT*b	C*a=c

C_cho.solve(c, x_cho);


// Singular Value Decomposition SVD
SVD svd(A);
VecDoub_O x_svd(3);

svd.solve(b,x_svd,svd.eps);     // Eps er threshold, med maskin præsision

// Calculating Standard Deviation
// sigma^2 = sum((V_ji/W_i)^2)
VecDoub std_dev(3);
Doub temp = 0;
for(int j=0;j<svd.n;j++){
    temp = 0;
    for(int i=0;i<svd.n;i++){
        temp += ((svd.v[j][i]/svd.w[i])*(svd.v[j][i]/svd.w[i]));
    }
    std_dev[j] = sqrt(temp);
}

//Det følgende virker - det ovenfor virker ikke???
//std_dev[0] = pow(svd.v[0][0]/svd.w[0],2) + pow(svd.v[0][1]/svd.w[1],2) + pow(svd.v[0][2]/svd.w[2],2);

//std_dev[1] = pow(svd.v[1][0]/svd.w[0],2) + pow(svd.v[1][1]/svd.w[1],2) + pow(svd.v[1][2]/svd.w[2],2);

//std_dev[2] = pow(svd.v[2][0]/svd.w[0],2) + pow(svd.v[2][1]/svd.w[1],2) + pow(svd.v[2][2]/svd.w[2],2);

//std_dev[0] = sqrt(std_dev[0]);

//std_dev[1] = sqrt(std_dev[1]);

//std_dev[2] = sqrt(std_dev[2]);


//util::print(A, "Matrix A");
//util::print(AT, "Transposed A");
//util::print(C, "C or alpha matrix AT*A");
//util::print(c, "c or beta vector AT*b");
//util::print(x, "LUdcmp result vector a, also called x");
//util::print(x_cho, "Cholesky result vector a, also called x");
//util::print(x_svd, "SVD result vector a, also called x");
util::print(std_dev, "Standard deviations");

// Checking result 
// calculate C*a = c
VecDoub c_lu = operator*(C,x);		//For LU dcmp
//util::print(c_lu, "c calculated from C*a - LU dcmp");
VecDoub c_cho = operator*(C, x_cho);	//For Cholesky
//util::print(c_cho, "c calculated from C*a - Cholesky");




}





void fili_data(){
VecDoub xFilip(82); VecDoub yFilip(82);
ifstream Filip("FilipData.dat");
for(int i = 0; i < 82; i++) {
	Filip >> yFilip[i];
	Filip >> xFilip[i];
}


Doub std_dev_fili = 0.003348;


// Normal Equation for solving problem
// (AT*A)*a = AT*b
// A design matrix,
// AT transpose of A
// a the vector to solve fore - typically named x
// b vector with result from equation A*x=b

// Creating Matrix A
MatDoub A(82,11);
for(int i=0;i<82;i++){
	for(int j=0;j<11;j++){
		A[i][j] = pow(xFilip[i],j)/std_dev_fili;
	}
}

// Creating vector b
VecDoub b(82);
for(int i=0;i<82;i++){
	b[i] = yFilip[i]/std_dev_fili;
}

// Transposing matrix A
MatDoub AT=util::Transpose(A);

// From normal equation
// AT*A = C
// Capital C is the resulting matrix
// In book capital C is alpha
//
// AT*b = c
// lowercase c is the resulting vector
// In book lowercase c is beta

// Creating Matrix C
MatDoub C= operator*(AT,A);

// Creating vector c
VecDoub c= operator*(AT,b);

// LU decomposition solve of linear equation problem
// This is supposed to fail because of numerical precission
LUdcmp lu(C);
VecDoub_O x(11);		//a in (AT*A)*a = AT*b	C*a=c

lu.solve(c, x);

// Cholesky solve of the linear equation problem
// This cannot compute since the A matrix is not positive definite - because of nummerical precission
//Cholesky C_cho(C);
//VecDoub_O x_cho(11);	//a in (AT*A)*a = AT*b	C*a=c

//C_cho.solve(c, x_cho);


// Singular Value Decomposition SVD
SVD svd(A);
VecDoub_O x_svd(11);

svd.solve(b,x_svd,svd.eps);

// Calculating Standard Deviation
// sigma^2 = sum((V_ji/W_i)^2)
VecDoub std_dev(11);
Doub temp = 0;
for(int j=0;j<svd.n;j++){
    for(int i=0;i<svd.n;i++){
       temp  += ((svd.v[j][i]/svd.w[i]) * (svd.v[j][i]/svd.w[i]));  
    }
    std_dev[j] = sqrt(temp);
}


//util::print(A, "Matrix A");
//util::print(AT, "Transposed A");
//util::print(C, "C or alpha matrix AT*A");
//util::print(c, "c or beta vector AT*b");
//util::print(x, "LUdcmp result vector a, also called x");
//util::print(x_cho, "Cholesky result vector a, also called x");
//util::print(x_svd, "SVD result vectora, also called x");
util::print(std_dev);


// Checking result 
// calculate C*a = c
//VecDoub c_lu = operator*(C,x);		//For LU dcmp
//util::print(c_lu, "c calculated from C*a - LU dcmp");
//VecDoub c_cho = operator*(C, x_cho);	//For Cholesky
//util::print(c_cho, "c calculated from C*a - Cholesky");		


}



int main() {
pont_data();
//fili_data();

return 0;
}
