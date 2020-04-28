#include "nr2.h"
#include "utilitires.h"
#include <iosteam>

using namespace std;


VecDoub initialGuess(int x0, int xn, int y0, int yn, int n){
    VecDout y(n);
    Doub a = (yn-y0)/(xn-x0);
    Doub b = y0-a*x0;
    for(int i = 0;i <n; i++){
            y[i] = a*i+b;
    }
    return y;
}

Doub F(dydx,y,x){

}


int main(){

    VecDoub y = initialGuess();
    
}
