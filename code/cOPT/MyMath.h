#ifndef _MYMATH_H_
#define _MYMATH_H_

#include "math.h"
#include <vector>

using namespace std;

extern double pINF;
#define MAX(x,y)	( ((x)>(y))?(x):(y) )
#define MIN(x,y)	( ((x)<(y))?(x):(y) )
#define ROUND(x)	( round( (x)*((double)1E12) )/((double)1E12) )
#define DELTA		(1E-3)

double lgamma(double x);
double lbeta(double n1,double n2);
double add_logged_values(double lx,double ly);
double sum_vector(vector<double>* x);

#endif


