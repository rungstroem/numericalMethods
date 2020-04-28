#include "nr3.h"
#include "utilities.h"
#include <iostream>
#include "rk4.h"
#include <fstream>
#include <math.h>
using namespace std;

Doub M = 5.68*pow(10.0,26.0);
Doub g = 4.98*pow(10.0,-10.0);
Doub m1 = 9.20*pow(10.0,18.0);
Doub m2 = 9.20*pow(10.0,18.0);

Doub r(Doub x, Doub y){
    return pow((pow(x,2.0) + pow(y,2.0)),0.5);
}

Doub R(Doub x1, Doub x2, Doub y1, Doub y2){
    return pow((pow((x1-x2),2.0) + pow((y1-y2),2.0)),0.5);
}


void derivs(Doub t, VecDoub_I &z, VecDoub_O &dzdt){
    dzdt[0] = z[4]; 
	dzdt[1] = z[5];
	dzdt[2] = z[6];
	dzdt[3] = z[7];
	dzdt[4] = -((M*g)/pow(r(z[0],z[1]),3)) * z[0] + ((m2*g)/pow(R(z[0],z[2],z[1],z[3]),3)) * (z[2]-z[0]); //x1
	dzdt[5] = -((M*g)/pow(r(z[0],z[1]),3)) * z[1] + ((m2*g)/pow(R(z[0],z[2],z[1],z[3]),3)) * (z[3]-z[1]); //y1
	dzdt[6] = -((M*g)/pow(r(z[2],z[3]),3)) * z[2] - ((m1*g)/pow(R(z[0],z[2],z[1],z[3]),3)) * (z[2]-z[0]); //x2
	dzdt[7] = -((M*g)/pow(r(z[2],z[3]),3)) * z[3] - ((m1*g)/pow(R(z[0],z[2],z[1],z[3]),3)) * (z[3]-z[1]); //y2
}

void convert(VecDoub &in, VecDoub &out){
    out[0] = pow( (pow(in[0],2.0) + pow(in[1],2.0)),0.5);
    out[1] = pow( (pow(in[2],2.0) + pow(in[3],2.0)),0.5);
    out[2] = (((atan2(in[1],in[0]))*(180.0/M_PI) - (atan2(in[3],in[2]))*(180.0/M_PI)));
}

VecDoub K(VecDoub results){
	//Order
    VecDoub K(results.size());
    K[0] = 0, K[1] = 0, K[2] = 0;
    Doub Ah1, Ah2, Ah3, ak2;

    for(int i = 0;i<results.size()-2;i++){
        Ah1 = results[i]; 
        Ah2 = results[i+1]; 
        Ah3 = results[i+2];
    
        ak2 = (Ah1-Ah2)/(Ah2-Ah3);
        K[i+2] = log(ak2)/log(2);
    }    
    return K;
}

VecDoub R_K(VecDoub results){
	//Order after Richardssons
    VecDoub K(results.size());
    K[0] = 0, K[1] = 0, K[2] = 0;
    Doub Ah1, Ah2, Ah3, ak2;

    for(int i = 1;i<results.size()-2;i++){
        Ah1 = results[i]; 
        Ah2 = results[i+1]; 
        Ah3 = results[i+2];
    
        ak2 = (Ah1-Ah2)/(Ah2-Ah3);
        K[i+2] = log(ak2)/log(2);
    }    
    return K;
}



VecDoub R_extrapolation(VecDoub results, int K, Doub alpha){
	//Richardssons extrapolation for optimizing results
    VecDoub Ri(results.size());
    Ri[0] = 0;
    Doub Ah1, Ah2;

    for(int i = 0;i<results.size()-1;i++){
        Ah1 = results[i];
        Ah2 = results[i+1];

        Ri[i+1] = (pow(alpha,K)*Ah2-Ah1)/(pow(alpha,2)-1);
    }
    return Ri;
    
}

VecDoub err_estimation(VecDoub results, int K, Doub alpha){
	//Richardssons error estimate
   VecDoub err(results.size());
   err[0] = 0;
   Doub Ah1, Ah2;

   for(int i = 0;i<results.size()-1;i++){
        Ah1 = results[i];
        Ah2 = results[i+1];

        err[i+1] = ((Ah1-Ah2)/(pow(alpha,K)-1));
   }
   return err;
}

VecDoub R_err_estimation(VecDoub results, int K, Doub alpha){
    VecDoub R_err(results.size());
    R_err[0] = 0, R_err[1] = 0;
    Doub Ah1, Ah2;

    for(int i = 1;i<results.size()-1;i++){
        Ah1 = results[i];
        Ah2 = results[i+1];

        R_err[i+1] = ((Ah1-Ah2)/(pow(alpha,K)-1));
    }
    return R_err;
}

int main(){
    
    Doub N = 500;
    Doub step = 1.0*pow(10.0,-3.0);
    Doub h = 0.0;

    VecDoub zinit(8);
    zinit[0] = 0.0; //x1
    zinit[1] = 152870.0;    //y1
    zinit[2] = 0.0; //x2
    zinit[3] = -153130.0;   //y2
    zinit[4] = -1360278.1;  //dx1
    zinit[5] = 0.0; //dy1
    zinit[6] = 1359122.8;   //dx2
    zinit[7] = 0;   //dy2
    VecDoub z(8);
    VecDoub_O dzdt(8);
    VecDoub_O zout(8);

    VecDoub polar(5);
	
    Doub t = 0.0;
    derivs(t, zinit, dzdt);
    z = zinit;
    t = 1.0;
    
    //File streaming
    ofstream dataset ("/home/runge/Desktop/yolo2.csv");
	
	VecDoub order(500000);
    VecDoub order2(500000);
    
    VecDoub result(10);
	VecDoub printttt(10);
    int i = 0;
    
    if(dataset.is_open()){
        for(double H = 1;H < 1024;) {
            
            for(; t < N; t+=1){
                derivs(t, z, dzdt);
                rk4(z, dzdt, t, h, zout, derivs);
                z = zout;
                convert(zout, polar);
                order[i] = polar[0];
                order2[i] = polar[1];
                i++;
            }
            
            //Polar coordinates
            convert(zout, polar);
            result[i] = polar[0];

		    //order[i] = polar[0];
            //order2[i] = polar[1];
		    
            H = H*2;
            h = 1/H;
            t = 1;
        }
    }
    
    printttt = K(result);
    util::print(printttt);
     
    if(dataset.is_open()){
        for(int j = 0;j<order.size();j++){
        dataset << result[j] << "," << order2[j] <<"\n";

        }
    }
    dataset.close();
    
    cout << "done\n";
    return 0;
}
