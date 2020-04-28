#include "nr3.h"
#include "utilities.h"
#include "ludcmp.h"
#include <iostream>
#include <cmath>

using namespace std;

Doub T1 = 1000.0;
Doub T2 = 500.0;
Doub eps1 = 0.8;
Doub eps2 = 0.6;
Doub sig = 1.712*pow(10,-9);
Doub d = 1.0;
Doub w = 1.0;       //Omega
Doub a = (-0.5)*w;  //Lower limit of integrals
Doub b = 0.5*w;     //Upper limit of integrals

Doub c1 = eps1*sig*pow(T1,4);
Doub c2 = eps2*sig*pow(T2,4);
Doub k1 = -(1-eps1);
Doub k2 = -(1-eps2);

Doub h = 0;

Doub f(Doub x, Doub y){
    return ( (1.0/2.0) * (1/pow( ( pow(d,2)+pow( (x-y),2 ) ),(3.0/2.0) ) ) );
}

template<class T>
VecDoub exercise2(int n, T &func){
    h = (b-a)/n;
    MatDoub A(2*(n+1), 2*(n+1));
    VecDoub B(2*(n+1));

    for(int i = 0 ; i<2*(n+1) ; i++){
        if( i<(n+1) ){
            B[i] = c1;
        }
        if( i>(n) ){
            B[i] = c2;
        }
    }

    for(int i = 0 ; i<2*(n+1) ; i++){
        for(int j =0 ; j<2*(n+1) ; j++){
            if(i==j){
                A[i][j] = 1;
            }
            else{
                A[i][j] = 0;   
            }
        }
    }

    int l = 0, g = 0;
    for(int i = 0 ; i<2*(n+1) ; i++){
        for(int j = 0; j<2*(n+1) ; j++){

            if( i<(n+1) && j>(n) ){
                A[i][j] = k1*h*func(a+i*h, a+l*h);
                l++;
            }
            if( i>n && j<(n+1) ){
                A[i][j] = k2*h*func(a+j*h, a+(i-(n+1))*h);
            }
        }
        l=0;
    }

    for(int i = 0;i<2*(n+1); i++){
        for(int j = 0;j<2*(n+1); j++){
            if( i<(n+1) && (j == n+1 || j == 2*n+1 )){
                A[i][j] *= 0.5;
            }
            if( i>n && (j == 0 || j == n)){
                A[i][j] *= 0.5;
            }
        }
    }

    //LU-decompostion
    VecDoub_O x_output(2*(n+1));
    LUdcmp lu(A);
    lu.solve(B,x_output);

    VecDoub test(2*(n+1));
    test = operator*(A, x_output);
    //util::print(x_output);
    return x_output;
    
}


//Task 2 determine Q1 & Q2
template<class T>
Doub Q1_func(VecDoub x_vec, T &func, int n){
	Doub hQ1 = (b-a)/n;
    Doub Q1 = 0;
    Doub temp = 0;
    Doub temp2 = 0;
    Doub temp3 = 0;
    Doub temp4 = 0;

    for(int i = 0; i<n+1; i++){
        for(int j = 0; j<n+1; j++){
            temp = h*func(a+i*h,a+j*h)*x_vec[j+n+1];
            if(j == 0 || j == n){
                temp2 += temp*0.5;
            }else{
                temp2 += temp;
            }
            temp = 0;
        }

        temp3 = h*(x_vec[i]-temp2);
        if(i==0 || i == n){
            temp4 += temp3*0.5;
        }else{
            temp4 += temp3;
        }
        
        temp3 = 0;
        temp2 = 0;

    }
    Q1 = temp4;
    
    /*
    VecDoub I1(n+1);
    for(int i = 0;i<n+1;i++){
        I1[i] = 0;
        for(int j = 0;j<n+1;j++){
            if(j == 0 || j == n){
                I1[i] += hQ1*0.5*func(a+hQ1*i, a+hQ1*j)*x_vec[j+n+1];
            }else{
                I1[i] += hQ1*func(a+hQ1*i,a+hQ1*j)*x_vec[j+n+1];
            }
        }
    }

    for(int i = 0;i<n+1;i++){
        if(i == 0 || i == n){
            Q1 += 0.5*hQ1*(x_vec[i]-I1[i]);
        }else{
            Q1 += hQ1*(x_vec[i]-I1[i]);
        }
    }
    */
    return Q1;
}

template<class T>
Doub Q2_func(VecDoub x_vec, T &func, int n){
	Doub Q2 = 0;
	Doub temp = 0;
    Doub temp2 = 0;
    Doub temp3 = 0;
    Doub temp4 = 0;

    for(int i = 0; i<n+1; i++){
        for(int j = 0; j<n+1; j++){
            temp = h*func(a+j*h,a+i*h)*x_vec[j];
            if(j == 0 || j == n){
                temp2 += temp*0.5;
            }else{
                temp2 += temp;
            }
            temp = 0;
        }
        temp3 = h*(x_vec[i+(n+1)]-temp2);
        if(i==0 || i==n){
            temp3 = temp3*0.5;
        }
        temp4 += temp3;
        temp3 = 0;
        temp2 = 0;

    }
    Q2 = temp4;
    return Q2;
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
    //const int n = 8;
    //VecDoub x(2*(n+1));
    //x = exercise2(n,f);
    //cout << Q1_func(x,f,n);


    const int n = 256;  //change iterations here!!
    VecDoub x(2*(n+1));

    VecDoub Q1_res(7);      //Change values here - if 4 iterations then value should be 1
    VecDoub Q2_res(7);      //If 8 ite. value should be 2, if 16 ite. value should be 3
    VecDoub Ri_Results1(7); //If 32 ite. value should be 4, if 64 ite. value should be 5
    VecDoub Ri_Results2(7); //If 128 ite. value should be 6, if 256 ite. value should be 7

    int j = 0;
    for(int i = 4;i<=n;){

        x = exercise2(i,f);
	    Q1_res[j] = Q1_func(x, f, i);
	    Q2_res[j] = Q2_func(x, f, i);
        i = i*2;
        j++;
    }

    /*
     * Hvad Virker
     * Exercise2 - dvs u og v
     * Q1
     * Q2
     * K af Q1
     * K af Q2
     * Richard extrapolation Q1
     * Richard extrapolation Q2
     * Richard ext. K1
     * Richard ext. K2
     *
     */
    cout << "u(x) & v(y) vector\n";
    util::print(x, "x-vector");
    cout << "\n";
    cout << "\n";
    cout << "\n";
    
    cout << "results for Q1\n";
    util::print(Q1_res, "Q1");
    cout << "\n";
    cout << "\n";
    util::print(K(Q1_res), "K for Q1");
    cout << "\n";
    cout << "\n";
    util::print(err_estimation(Q1_res,2,2), "Err est. af Q1");
    cout << "\n";
    cout << "\n";
    Ri_Results1 = R_extrapolation(Q1_res,2,2);
    util::print(Ri_Results1, "Ric. ext. Q1");
    cout << "\n";
    cout << "\n";
    util::print(R_K(Ri_Results1), "K for ric. ext. Q1");
    cout << "\n";
    cout << "\n";
    util::print(R_err_estimation(Ri_Results1,4,2), "Err est. af ric. ext. Q1");
    cout << "\n";
    cout << "\n";
    cout << "\n";
    cout << "results for Q2\n";
    util::print(Q2_res, "Q2");
    cout << "\n";
    cout << "\n";
    util::print(K(Q2_res), "K for Q2");
    cout << "\n";
    cout << "\n";
    util::print(err_estimation(Q2_res,2,2),"Err est. af Q2");
    cout << "\n";
    cout << "\n";
    Ri_Results2 = R_extrapolation(Q2_res,2,2);
    util::print(Ri_Results2, "Ric. ext. Q2");
    cout << "\n";
    cout << "\n";
    util::print(R_K(Ri_Results2), "K for ric. ext. Q2");
    cout << "\n";
    cout << "\n";
    util::print(R_err_estimation(Q1_res, 2,2), "Err est. af ric. ext. Q2");  

}
