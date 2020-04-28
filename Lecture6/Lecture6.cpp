#include "nr3.h"
#include "utilities.h"
#include "roots.h"
#include <math.h>
#include <iostream>

Doub pi = 3.1415926535;
Doub x1 = 0, xh = pi/2;

Doub f(const Doub x){
    return ( x-cos(x) );
}

//Use functor for defining the function when using Newton-method...

void bisection(){
    Doub C = 0.5; //Convergence constant
    Doub xacc = 1*pow(10,-8);
    Doub result = 0;
    MatDoub bisTable(50,6);

    result = rtbis(f, x1, xh, xacc);
    cout << "Bisection result\n";
    cout << result << "\n";
    //Result Checking
    cout << f(result) << "\n";
    
    //For Table
    bisTable = rtbisM(f,x1,xh,xacc);
    util::print(bisTable, "Bisection table");
}

void secant(){
    Doub C; //Convergence constant
    Doub result = 0, xacc = 1*pow(10,-16);
    MatDoub secTable(30,6);

    result = rtsec(f,x1,xh,xacc);
    cout << "Secant result\n";
    cout << result << "\n";

    //Result Checking
    cout << f(result) << "\n";

    //For Table
    secTable = rtsecM(f,x1,xh,xacc);
    util::print(secTable, "Secant table");
}

void riddle(){
    Doub result = 0, xacc = 1*pow(10,-16);
    MatDoub ridTable(60,4);

    result = zriddr(f,x1,xh,xacc);
    cout << "Riddle result\n";
    cout << result << "\n";

    //For Table
    ridTable = zriddrM(f,x1,xh,xacc);
    util::print(ridTable, "Riddle table");
}

int main (){
    //bisection();
    secant();
    //riddle();
}
