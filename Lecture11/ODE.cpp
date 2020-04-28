#include "nr3.h"
#include "utilities.h"
#include <iostream>

using namespace std;

VecDoub f(Doub x, Doub y){
    VecDoub result(2);

    Doub dxdt = x*y;
    Doub dydt = -(pow(x,2));

    result[0] = dxdt, result[1] = dydt;
    return result;

}

VecDoub K(VecDoub results){
	VecDoub K(results.size());
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


template<class T>
VecDoub euler(T &func, Doub t, Doub h, VecDoub yn){
    //Called Runge-Kutta 1'st order
    VecDoub yn1(2), temp(2);
    temp = func(yn[0],yn[1]);

    yn1[0] = yn[0] + h*temp[0];
    yn1[1] = yn[1] + h*temp[1];

    return yn1;
}
/*
template<class T>
VecDoub midpoint(T $func, Doub t, Doub h, VecDoub yn){
    //Called Runge-Kutta 2'nd order
    VecDoub yn1(2), k1(2), k2(2);
    k1 = func(yn[0], yn[1]);
    k1[0] = k1[0]*h;
    k1[1] = k1[1]*h;
    
    k2 = func(yn[0]+0.5*h,yn[1]+0,5*k1[0]);

    return k2;

}
*/

int main(){
    int n = 10;
    Doub h = 1.0/5.0;
    Doub x0 = 1.0, y0 = 1.0;

    VecDoub yn(2), yn1(2); 
	VecDoub results1(7);
	VecDoub results2(7);

    yn[0] = x0, yn[1] = y0;
	int j = 0;
	for(int i = 5;i<640+1;){
		for(int t = 0;t<n+1;t++){
			yn1 = euler(f,t,pow(i,-1),yn);
			yn = yn1;
		}
		results1[j] = yn1[0];
		results2[j] = yn1[1];
		j++;
		i = i*2;
	}
	util::print(results1);
	util::print(K(results1));
	VecDoub RiResults(7);
	RiResults = R_extrapolation(results1, 1,2);
	util::print(RiResults);
	cout << "\n";
	cout << "\n";
	util::print(results2);
	util::print(K(results2));
	RiResults = R_extrapolation(results2, 1, 2);
	util::print(RiResults);

}
