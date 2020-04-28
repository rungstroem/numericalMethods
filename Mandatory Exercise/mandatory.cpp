#include "nr3.h"
#include "utilities.h"
#include "svd.h"
#include <fstream>
#include <math.h>

using namespace std;

void dataset1(){

//////////////     Task 2!!!!! /////////////////
cout << "Task2\n";
//Creating 4 vectors for data import
VecDoub Theta1(500), Theta2(500), x(500), y(500); 

//Importing data
ifstream data1("d1");
for(int i = 0; i<500; i++){
    data1 >> Theta1[i];
    data1 >> Theta2[i];
    data1 >> x[i];
    data1 >> y[i];
}

int row = 1000, col = 4;
MatDoub A(1000,4);
VecDoub z(1000);
int j = 0;

for(int i=0; i<row; i+=2){
    //Constructing the A matrix
    A[i][0] = 1;
    A[i][1] = 0;
    A[i][2] = cos(Theta1[j]);
    A[i][3] = cos(Theta1[j]+Theta2[j]);

    A[i+1][0] = 0;
    A[i+1][1] = 1;
    A[i+1][2] = sin(Theta1[j]);
    A[i+1][3] = sin(Theta1[j]+Theta2[j]);
    
    //Constructing the z vector
    z[i] = x[j];
    z[i+1] = y[j];
    j++;
}
//util::print(A, "Matrix A");
//util::print(z, "Vector z");

//Doing the SVD decomposition
SVD svd(A);

//Printing U, W and V
util::print(svd.u , "U matrix");
cout << "\n";
util::print(svd.w , "W matrix");
cout << "\n";
util::print(svd.v , "V matrix");
cout << "\n";

// Da W matricens diagonal har værdier højere end 0 (positiv) kan vi sige at
// der ikke er lineær dependency i A matricen



////////////   Task 3!!! ////////
cout << "Task 3\n";
//Doing the SVD calculation
//Estimate q
VecDoub_O q(4);
svd.solve(z,q,svd.eps);

//Printing the result from SVD calculation
util::print(q, "Parameter estimation - the SVD result of q-vector");
cout << "\n";
//The residual error - ||Aq-z||
VecDoub_O res_err(1000);
res_err = operator*(A,q);
for(int i=0;i<1000;i+=1){
    res_err[i] = abs(res_err[i]-z[i]);
}

//Printing the residual error
util::print(res_err, "The Residual errors");
cout << "\n";


////////////   Task 4!!! ////////
cout << "Task4\n";
//The standard measurement error is 1mm - 0.1cm
double std_dev = 0.1;

//Creating a new A matrix with the standard deviation is incorporated
MatDoub A_kor(1000,4);
VecDoub z_kor(1000);
j = 0;
for(int i=0; i<row; i+=2){
    //Constructing the A matrix
    A_kor[i][0] = 1/std_dev;
    A_kor[i][1] = 0/std_dev;
    A_kor[i][2] = cos(Theta1[j])/std_dev;
    A_kor[i][3] = cos(Theta1[j]+Theta2[j])/std_dev;
    A_kor[i+1][0] = 0/std_dev;
    A_kor[i+1][1] = 1/std_dev;
    A_kor[i+1][2] = sin(Theta1[j])/std_dev;
    A_kor[i+1][3] = sin(Theta1[j]+Theta2[j])/std_dev;
    //Constructing the z vector
    z_kor[i] = x[j]/std_dev;
    z_kor[i+1] = y[j]/std_dev;
    j++;
}

SVD svd_kor(A_kor);
VecDoub_O q_kor(4);
svd.solve(z_kor, q_kor, svd_kor.eps);

util::print(q_kor, "Korrected q vector - with 0.1 cm std.dev");
cout << "\n";

// Calculating Standard Deviation
// sigma^2 = sum((V_ji/W_i)^2)
VecDoub std_dev_out(4);
Doub john = 0;
for(int j=0;j<4;j++){
    john = 0;
    for(int i=0;i<svd_kor.n;i++){
        john += ((svd_kor.v[j][i]/svd_kor.w[i])*(svd_kor.v[j][i]/svd_kor.w[i]));
    }
    std_dev_out[j] = sqrt(john);
}

//Estimated resulting error
util::print(std_dev_out, "Estimated resulting error - with EQ 15.4.19");
cout << "\n";

//////// Task 5!!! ////////
cout << "Task5\n";
//Calculate the mean for residuals
double x_bar = 0;
for(int i = 0;i<1000;i++){
    x_bar = x_bar + res_err[i];
}
x_bar = x_bar/1000.0;
cout << "\n";
cout << "Mean dataset 1 : " << x_bar << "\n";

//Calculate the Variance for residuals
double s2 = 0;
double temp2 = 0;
for(int i=0;i<1000;i++){
    temp2 = temp2 + pow((res_err[i]-x_bar),2);
}

s2 = temp2/(1000-1);
cout << "Variance dataset 1 : " << s2 << "\n";
}

void dataset2(){

//////////////     Task 2!!!!! /////////////////

//Creating 4 vectors for data import
VecDoub Theta1(500), Theta2(500), x(500), y(500); 

//Importing data
ifstream data1("d2");
for(int i = 0; i<500; i++){
    data1 >> Theta1[i];
    data1 >> Theta2[i];
    data1 >> x[i];
    data1 >> y[i];
}

int row = 1000, col = 4;
MatDoub A(1000,4);
VecDoub z(1000);
int j = 0;

for(int i=0; i<row; i+=2){
    //Constructing the A matrix
    A[i][0] = 1;
    A[i][1] = 0;
    A[i][2] = cos(Theta1[j]);
    A[i][3] = cos(Theta1[j]+Theta2[j]);

    A[i+1][0] = 0;
    A[i+1][1] = 1;
    A[i+1][2] = sin(Theta1[j]);
    A[i+1][3] = sin(Theta1[j]+Theta2[j]);
    
    //Constructing the z vector
    z[i] = x[j];
    z[i+1] = y[j];
    j++;
}
//util::print(A, "Matrix A");
//util::print(z, "Vector z");

//Doing the SVD decomposition
SVD svd(A);

//Printing U, W and V
//util::print(svd.u , "U matrix");
//util::print(svd.w , "W matrix");
//util::print(svd.v , "V matrix");

// Da W matricens diagonal har værdier højere end 0 (positiv) kan vi sige at
// der ikke er lineær dependency i A matricen



////////////   Task 3!!! ////////

//Doing the SVD calculation
//Estimate q
VecDoub_O q(4);
svd.solve(z,q,svd.eps);

//Printing the result from SVD calculation
//util::print(q, "The SVD result of q-vector");

//The residual error
VecDoub_O res_err(1000);
res_err = operator*(A,q);
for(int i=0;i<1000;i+=1){
    res_err[i] = abs(res_err[i]-z[i]);
}

//Printing the residual error
//util::print(res_err, "The Residual errors");



////////////   Task 4!!! ////////

//The standard measurement error is 1mm - 0.1cm
double std_dev = 0.1;

//Creating a new A matrix with the standard deviation is incorporated
MatDoub A_kor(1000,4);
VecDoub z_kor(1000);
j = 0;
for(int i=0; i<row; i+=2){
    //Constructing the A matrix
    A_kor[i][0] = 1/std_dev;
    A_kor[i][1] = 0/std_dev;
    A_kor[i][2] = cos(Theta1[j])/std_dev;
    A_kor[i][3] = cos(Theta1[j]+Theta2[j])/std_dev;
    A_kor[i+1][0] = 0/std_dev;
    A_kor[i+1][1] = 1/std_dev;
    A_kor[i+1][2] = sin(Theta1[j])/std_dev;
    A_kor[i+1][3] = sin(Theta1[j]+Theta2[j])/std_dev;
    //Constructing the z vector
    z_kor[i] = x[j]/std_dev;
    z_kor[i+1] = y[j]/std_dev;
    j++;
}

SVD svd_kor(A_kor);
VecDoub_O q_kor(4);
svd.solve(z_kor, q_kor, svd_kor.eps);

//util::print(q_kor, "korrected q vector");

// Calculating Standard Deviation
// sigma^2 = sum((V_ji/W_i)^2)
VecDoub std_dev_out(4);
Doub john = 0;
for(int j=0;j<4;j++){
    john = 0;
    for(int i=0;i<svd_kor.n;i++){
        john += ((svd_kor.v[j][i]/svd_kor.w[i])*(svd_kor.v[j][i]/svd_kor.w[i]));
    }
    std_dev_out[j] = sqrt(john);
}

//Estimated resulting error
util::print(std_dev_out, "Standard Deviation korrected results");


//////// Task 5!!! ////////

//Calculate the mean for residuals
double x_bar = 0;
for(int i = 0;i<1000;i++){
    ;x_bar = x_bar + res_err[i];
}
x_bar = x_bar/1000.0;
cout << "\n";
cout << "Mean dataset 2 : " << x_bar << "\n";

//Calculate the Variance for residuals
double s2 = 0;
double temp2 = 0;
for(int i=0;i<1000;i++){
    temp2 = temp2 + pow((res_err[i]-x_bar),2);
}

s2 = temp2/(1000-1);
cout << "Variance dataset 2 : " << s2 << "\n";
}


int main(){
    dataset1();
    dataset2();
    return 0;

}
