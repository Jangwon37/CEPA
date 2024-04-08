#ifndef _WKM_H_
#define _WKM_H_

#include <vector>
#include "CenData.h"
#include "Density.h"
#include "Region.h"

using namespace std;

class WKM
{
public:
	WKM(vector<cendata_t>* data);
	WKM(vector<cendata_t>* data,vector<double>* weights);
	~WKM() {};

	Density* getDensity() { return &mWKMDensity; }
	double getSurvProb(double t);
	double getCondSurvProb(double t,double t0);
	double getProb(double t1, double t2) { return getSurvProb(t1)-getSurvProb(t2); }

private:
	void calculateWeightedKaplanMeier(vector<cendata_t>* data,vector<double>* weights);
	Density mWKMDensity;
	vector<double> mTime;
	vector<double> mProb;

};

#endif


