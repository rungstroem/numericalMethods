#include "nr3.h"
#include "utilities.h"
#include <iostream>
#include "ludcmp.h"
#include "qrdcmp.h"
#include "roots_multidim.h"

using namespace std;

bool CheckFlag;
VecDoub_IO vec_x(8);    //x-vector
//VecDoub_IO x_vec(8);

VecDoub vecfunc(VecDoub_I x_vec){

	//Konstanter
	Doub v = 120;
	Doub K = 2.5;
	Doub w = 4.0;
	Doub alpha = 2*pow(10,-7);
	Doub d = 30;
	Doub n = 0.1;
	
	//Variables
	Doub Lo = x_vec[0];
	Doub L = x_vec[1];
	Doub p = x_vec[2];
	Doub x = x_vec[3];
	Doub theta = x_vec[4];
	Doub phi = x_vec[5];
	Doub a = x_vec[6];
	Doub H = x_vec[7];
    
	//YoloSwag
	VecDoub results(8);

	results[0] = L/(1+alpha*H);			//Lo
	results[1] = 2*a*sinh(x/a);			//L
	results[2] = n-K*sin(theta);			//p
	results[3] = (d-2*K*cos(theta))/2;		//x
	results[4] = atan((1+(v/(w*Lo)))*tan(phi));	//theta
	results[5] = atan(sinh(x/a));			//phi
	results[6] = p/(cosh(x/a)-1);			//a
	results[7] = (w*Lo)/(2*sin(phi));		//H
	
	//Return vector of functions...
	return results;
}

int main(){
//Initial variables
Doub Lo_init = 29;
Doub L_init = 25;
Doub p_init = 2;
Doub x_init = 14;
Doub theta_init = 1.92;
Doub phi_init = 0.35;
Doub H_init = 5;
Doub a_init = 40;

//Setting initial variables
vec_x[0] = Lo_init;
vec_x[1] = L_init;
vec_x[2] = p_init;
vec_x[3] = x_init;
vec_x[4] = theta_init;
vec_x[5] = phi_init;
vec_x[6] = a_init;
vec_x[7] = H_init;


//Call of Newton function
newt(vec_x,CheckFlag,vecfunc);

}
