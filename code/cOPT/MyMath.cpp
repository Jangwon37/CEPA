#include <math.h>
#include "MyMath.h"

// lgamma.cpp -- log gamma function of real argument.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
//  Returns log(gamma) of real argument.
//  NOTE: Returns 1e308 if argument is 0 or negative.
double lgamma(double x) 
{
	double x0,x2,xp,gl,gl0;
	int n=0,k;
	static double a[] = {
		 8.333333333333333e-02,
       		-2.777777777777778e-03,
		 7.936507936507937e-04,
		-5.952380952380952e-04,
		 8.417508417508418e-04,
		-1.917526917526918e-03,
		 6.410256410256410e-03,
		-2.955065359477124e-02,
		 1.796443723688307e-01,
		-1.39243221690590};

	x0 = x;
	//if (x <= 0.0) return 1e308;
	if (x <= 0.0) return INFINITY;
	else if ((x == 1.0) || (x == 2.0)) return 0.0;
	else if (x <= 7.0) {
		n = (int)(7-x);
		x0 = x+n;
	}
	x2 = 1.0/(x0*x0);
	xp = 2.0*M_PI;
	gl0 = a[9];
	for (k=8;k>=0;k--) {
		gl0 = gl0*x2 + a[k];
	}
	gl = gl0/x0+0.5*log(xp)+(x0-0.5)*log(x0)-x0;
	if (x <= 7.0) {
		for (k=1;k<=n;k++) {
			gl -= log(x0-1.0);
			x0 -= 1.0;
		}
	}
    
	return gl;
}

// log beta function
double lbeta(double x,double y) 
{
	return( lgamma(x)+lgamma(y)-lgamma(x+y) );
}

// return log(x+y) by given log(x) and log(y)
double add_logged_values(double lx,double ly)
{
	double common, left;
	if( isinf(lx) ) return ly;
	if( isinf(ly) ) return lx;
	common = MAX(lx,ly);
	left = MIN(lx,ly) - common;
	return( common + log(1+exp(left)) );
}

// sum a vector
double sum_vector(vector<double>* x)
{
	double s=0.0;
	for( unsigned int i=0 ; i<x->size() ; i++ ) s += x->at(i);
        return s;
}


