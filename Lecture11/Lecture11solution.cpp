#include <iostream>
#include <nr3.h>
#include <stepper.h>
#include <stepperdopr853.h>
#include <stepperdopr5.h>
#include <odeint.h>
#include <shoot.h>

#include "quadrature.h"
#include "qrdcmp.h"
#include "ludcmp.h"
#include "svd.h"
#include "roots_multidim.h"
#include "psplot.h"
using namespace std;


struct rhs_van {
	Doub eps;
	rhs_van(Doub epss) : eps(epss) {}
	void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx)
	{
		dydx[0] = y[1];
		dydx[1] =-cos(y[0])*sin(y[1]);
	}
};

struct load {
	load() {}
	VecDoub operator() (const Doub x, VecDoub_I &y)
	{
		VecDoub ystart(2);
		ystart[0] = 0;
		ystart[1] = y[0];
		return ystart;
	}
};

struct score {
	Doub t;
	score(Doub target) :t(target) {}
	VecDoub operator() (const Doub x, VecDoub_I &y)
	{
		VecDoub error(1);
		error[0] =  y[0] - t;
		 return error;
	}
};


int main() {

	Int nvar = 2;
	Doub x1 = 0.0, x2 = 10.0;
	load l = load();
	rhs_van dd(10e-3);
	score s(3);
	VecDoub v(1);

	bool check = false;

	cout << endl << "Shooting to find initial conditions" << endl;
	Shoot<load,rhs_van,score> shoot(nvar,x1,x2,l,dd,s);
	newt(v,check,shoot);



	if(!check)
	{
		cout << v[0] << endl;
	}

	cout << endl << "Find solution using stepperdorph5" << endl;
	const Int nnvar = 2;
	const Doub 	atol=1.0e-14,
				rtol=atol,
				h1=(x2-x1)/100.0,
				hmin=0.0;

	VecDoub ystart(nnvar);
	ystart[0]=0;
	ystart[1]=v[0];

	Output outDorph5(20);
	rhs_van d(10e-14);
	try {
		Odeint<StepperDopr5<rhs_van> > odeDopr5(ystart,x1,x2,atol,rtol,h1,hmin,outDorph5,d);
		odeDopr5.integrate();
	} catch (NRerror e) {cout << e.message << endl;}

	cout << "outDorph5.count: " << outDorph5.count << endl;
	for(Int i =0; i<outDorph5.count;i++)
	{
		cout << setprecision(8) << outDorph5.xsave[i] << "\t" << outDorph5.ysave[0][i] << "\t" << outDorph5.ysave[1][i] << endl;
	}

	cout << "End" << endl; // prints the world is beautiful

	PSpage pg("myplot.ps");
	PSplot plot1(pg,100.,500.,500.,800.);
	//PSplot plot2(pg,100.,500.,100.,400.);
	plot1.setlimits(x1,x2,0.,5.);
	plot1.frame();
	plot1.autoscales();
	plot1.xlabel("x");
	plot1.ylabel("y");
	//plot1.lineplot(xlag,ylag);

	for(int i = 0; i< 21; i++)
	{
		plot1.pointsymbol(outDorph5.xsave[i],outDorph5.ysave[0][i],108,8);
	}

	for(int i = 0; i< 21; i++)
	{
		plot1.pointsymbol(outDorph5.xsave[i],outDorph5.ysave[1][i],100,8);
	}

	pg.close();
	pg.display("\"ghostscript\" ");


	return 0;
}
