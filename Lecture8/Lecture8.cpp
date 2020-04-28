#include "nr3.h"
#include "utilities.h"
#include <iostream>
#include <cmath>

Doub e = M_E; 

Doub f1(Doub x){
    return cos(pow(x,2))*pow(e,-x);
}

Doub f2(Doub x){
    return pow(x,0.5)*cos(pow(x,2))*pow(e,-x);
}

Doub f3(Doub x){
    return 1000*pow(e,(-1/x))*pow(e,(-1/(1-x)));
}

Doub f4(Doub x){
    return (1/pow(x,0.5))*cos(pow(x,2))*pow(e,-x);
}

Doub f5(Doub x){
    return (pow(e,x));
}

template<class T>
Doub Midpoint(Doub x0, Doub xn, int n, T &func){
    Doub h = (xn-x0)/n;
    Doub Si = 0;
    Doub Sf = 0;
    for(int i = 0;i<n;i++){ 
       Si = h*func(x0+i*h);
       Sf += Si;
       Si = 0;
    }

    return Sf;
}

template<class T>
VecDoub Midpoint2(Doub x0, Doub xn, int n, T &func){
    Doub h = (xn-x0)/n;
    Doub Si = 0;
    Doub Sf = 0;

    VecDoub results(n);

    for(int i = 0;i<n;i++){ 
       Si = h*func(x0+i*h);
       Sf += Si;
       results[i] = Sf;
       Si = 0;
    }

    return results;
}

template<class T>
Doub Trapezoidal(Doub x0, Doub xn, int n, T &func){
    Doub h = (xn-x0)/n;
    Doub Sf = h*0.5*(func(x0)+func(xn));
    Doub Si = 0;

    for(int i = 1; i < n; i++){
         
        Si = h*func(x0+i*h);
        Sf += Si;
        
    }

    return Sf;
}

template<class T>
VecDoub Trapezoidal2(Doub x0, Doub xn, int n, T &func){
    Doub h = (xn-x0)/n;
    Doub Sf = h*0.5*(func(x0)+func(xn));
    Doub Si = 0;
    VecDoub results(n);

    for(int i = 1; i < n; i++){

        Si = h*func(x0+i*h);
        Sf += Si;
        results[i] = Sf;
        
    }

    return results;
}

template<class T>
Doub Simpsons(Doub x0, Doub xn, int n, T &func){
    Doub h = (xn-x0)/n;
    Doub Si = 0;
    Doub Sf = 0;
    Doub S = 0;
    for(int i = 1; i<n;i++){
        if(i%2==0){
            Si += 2*func(x0+i*h);
        }
        else{
            Sf += 4*func(x0+i*h);
        }
    }
    
    S = (1.0/3.0)*h*(func(x0)+func(xn)+Si+Sf);
    
    return S;
} 

VecDoub K(VecDoub results){
    VecDoub K(results.size());
    K[0] = 0, K[1] = 0;
    Doub Ah1, Ah2, Ah3, ak2;

    for(int i = 0;i<results.size()-2;i++){
        Ah1 = results[i]; 
        Ah2 = results[i+1]; 
        Ah3 = results[i+2];
    
        ak2 = (Ah1-Ah2)/(Ah2-Ah3);
        K[i+2] = log2(ak2);
    }    
    return K;
}

VecDoub R_extrapolation(VecDoub results, int K, Doub alpha){
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

VecDoub R_err_estimation(VecDoub results, int K, Doub alpha){
   VecDoub Ri_err(results.size());
   Ri_err[0] = 0;
   Doub Ah1, Ah2;

   for(int i = 0;i<results.size()-1;i++){
        Ah1 = results[i];
        Ah2 = results[i+1];

        Ri_err[i+1] = ((Ah1-Ah2)/(pow(alpha,K)-1));
   }
   return Ri_err;
}

int main(){
    int n = 32;
    VecDoub results(n/5);
    VecDoub Ri_Results(n/5);
    int j = 0;
    for(int i = 1;i<=n;){
        results[j] = Trapezoidal(0,4,i,f5);
        i = i*2;
        j++;
    }
    util::print(results);
    cout << "\n";
    cout << "\n";
    util::print(K(results));
    cout << "\n";
    cout << "\n";
    Ri_Results = R_extrapolation(results,2,2);
    util::print(Ri_Results);
    cout << "\n";
    cout << "\n";
    util::print(K(Ri_Results));
    cout << "\n";
    cout << "\n";
    util::print(R_err_estimation(results, 2,2));
    //cout << Midpoint(0,4,2,f5) << "\n";
    //cout << Trapezoidal(0,4,32,f5) << "\n";
    //cout << Simpsons(0,4,2,f5) << "\n";
}
