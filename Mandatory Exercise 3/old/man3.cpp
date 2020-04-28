#include "nr3.h"
#include "utilities.h"
#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

Doub M = 5.68*pow(10,26);
Doub g = 4.98*pow(10,-10);
Doub m1 = 9.2*pow(10,18);
Doub m2 = 9.2*pow(10,18);

//some values
void initial(VecDoub &y){
    //Start værdier
    // x & y
    y[0] = 0.0;
    y[1] = 152870;
    y[2] = 0.0;
    y[3] = -153130;

    //dydx
    y[4] = -1360278.1;
    y[5] = 0.0;
    y[6] = 1359122.8;
    y[7] = 0.0;
}

Doub r(Doub x, Doub y){
    return(pow(pow(x,2)+pow(y,2),0.5));
}

Doub r12(Doub x1, Doub x2, Doub y1, Doub y2){
    return(pow( pow(x1-x2,2)+pow(y1-y2,2), 0.5));
}
/*
//some shady functions
Doub Zx1(Doub x1){
    return x1;
}
Doub Zx2(Doub x2){
    return x2;
}
Doub Zy1(Doub y1){
    return y1;
}
Doub Zy2(Doub y2){
    return y2;
}
Doub dZx1(Doub x1, Doub x2, Doub y1, Doub y2, Doub M, Doub m2, Doub g){
    return (((M*g)/pow(r(x1,y1),3)) *x1+((m2*g)/pow(r12(x1,x2,y1,y2),3))*(x2-x1));
}
Doub dZx2(Doub x1, Doub x2, Doub y1, Doub y2, Doub M, Doub m1, Doub g){
    return (((M*g)/pow(r(x2,y2),3)) *x2-((m1*g)/pow(r12(x1,x2,y1,y2),3))*(x2-x1));
}
Doub dZy1(Doub x1, Doub x2, Doub y1, Doub y2, Doub M, Doub m2, Doub g){
    return (((M*g)/pow(r(x1,y1),3)) *y1+((m2*g)/pow(r12(x1,x2,y1,y2),3))*(y2-y1));
}
Doub dZy2(Doub x1, Doub x2, Doub y1, Doub y2, Doub M, Doub m1, Doub g){
    return (((M*g)/pow(r(x2,y2),3)) *y2-((m1*g)/pow(r12(x1,x2,y1,y2),3))*(y2-y1));
}

void derivs(Doub t, VecDoub_I &init, VecDoub_O &out){
    out[0] = Zx1(init[4]);
    out[2] = Zx2(init[5]);
    out[1] = Zy1(init[6]);
    out[3] = Zy2(init[7]);
    out[4] = dZx1(init[0], init[2], init[1], init[3], M, m2, g);
    out[6] = dZx2(init[0], init[2], init[1], init[3], M, m1, g);
    out[5] = dZy1(init[0], init[2], init[1], init[3], M, m2, g);
    out[7] = dZy2(init[0], init[2], init[1], init[3], M, m1, g);
}
*/
void derivs(Doub t, VecDoub_I &init, VecDoub_O &out){

    out[0] = init[4];
    out[1] = init[5];
    out[2] = init[6];
    out[3] = init[7];
   	out[4] = -((M*g)/pow(r(init[0],init[1]),3))*init[0] + ((m2*g)/pow(r12(init[0], init[1], init[2], init[3]),3))*(init[2]-init[0]);
	out[5] = -((M*g)/pow(r(init[0],init[1]),3))*init[1] + ((m2*g)/pow(r12(init[0], init[1], init[2], init[3]),3))*(init[3]-init[1]);
	out[6] = -((M*g)/pow(r(init[0],init[1]),3))*init[2] - ((m2*g)/pow(r12(init[0], init[1], init[2], init[3]),3))*(init[2]-init[0]);
	out[7] = -((M*g)/pow(r(init[0],init[1]),3))*init[3] - ((m2*g)/pow(r12(init[0], init[1], init[2], init[3]),3))*(init[3]-init[1]);
	

}

//All credit goes to me.
void rk4(VecDoub_I &y, VecDoub_I &dydx, const Doub x, const Doub h, VecDoub_O &yout, void derivs(const Doub, VecDoub_I &, VecDoub_O &)) {

        Int n=y.size();
        VecDoub dym(n),dyt(n),yt(n);
        Doub hh=h*0.5;
        Doub h6=h/6.0;
        Doub xh=x+hh;
        for (Int i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
        derivs(xh,yt,dyt);
        for (Int i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
        derivs(xh,yt,dym);
        for (Int i=0;i<n;i++) {
                yt[i]=y[i]+h*dym[i];
                dym[i] += dyt[i];
        }
        derivs(x+h,yt,dyt);
        for (Int i=0;i<n;i++)
                yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
}

MatDoub polar(VecDoub in){
    
    MatDoub out(2,2);
    out[0][0] = pow(pow(in[0],2)+pow(in[1],2),0.5);
    out[0][1] = atan2(in[1],in[0]);
    out[1][0] = pow(pow(in[2],2)+pow(in[3],2),0.5);
    out[1][1] = atan2(in[3],in[2]);

    return out;

}


int main(){
	Doub N = 500; //Antal date
	Doub h = 1*pow(10,-5);//(b-0.0)/N;
    
    VecDoub y(8);
    initial(y);
	VecDoub dydx(8);

	VecDoub_O yout(8);
    

    //Polar out
    VecDoub useless(2);
    MatDoub fuck(2,2);
    ofstream dataset ("/home/runge/Desktop/data.csv");
    ofstream dataset2 ("/home/runge/Desktop/data2.csv");
    ofstream dataset3 ("/home/runge/Desktop/data3.csv");
    if(dataset.is_open() && dataset2.is_open()){
	for(Doub t = 0.0; t<N; t+=h){
        derivs(t,y,dydx);
		rk4(y, dydx, t, h, yout, derivs);
	    y = yout;
        
        fuck = polar(yout);
        useless[0] = abs(fuck[0][0] - fuck[1][0]);
        useless[1] = abs(fuck[0][1] - fuck[1][1]);

        //util::print(polar(yout));
        dataset << t << "," << fuck[0][0] << "," << fuck[0][1] <<  "\n";
        dataset2 << t << "," << fuck[1][0] << "," << fuck[1][1] << "\n";
        dataset3 << t << "," << useless[1] << "\n";
	}
    }
    dataset.close();
    	
    return 0;
}
