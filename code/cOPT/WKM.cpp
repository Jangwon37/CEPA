#include <iostream>
#include <algorithm>
#include <math.h>
#include "WKM.h"
#include "MyMath.h"

using namespace std;

WKM::WKM(vector<cendata_t>* data) 
{
	vector<double> weights(data->size(),1.0);
	calculateWeightedKaplanMeier(data,&weights);
}

WKM::WKM(vector<cendata_t>* data,vector<double>* weights) 
{
	calculateWeightedKaplanMeier(data,weights);
}

// weighted kaplan meier
void WKM::calculateWeightedKaplanMeier(vector<cendata_t>* data,vector<double>* weights)
{
	// collected non-censored time points and sort
	vector<double> nct; bool has_zero = false; 
	//double max_time = -1.0; bool censored_max_time = false;
	for( unsigned int i=0 ; i<data->size() ; i++ ) {
		if( data->at(i).e == 1 ) {
			if( data->at(i).x == 0.0 ) has_zero = true;
			nct.push_back(data->at(i).x);
	//		if( data->at(i).x > max_time ) {
	//			max_time = data->at(i).x;
	//			if( data->at(i).e == 1 ) censored_max_time = false;
	//			else censored_max_time = true;
	//		}
		}
	}
	if( !has_zero ) nct.push_back(0.0);
	sort(nct.begin(),nct.end());

	// estimate a survival function
	vector<double> surv_prob(nct.size(),1);
	for( unsigned int i=1 ; i<nct.size() ; i++ ) {
		double r = 0.0, d = 0.0;
		for( unsigned int j=0 ; j<data->size() ; j++ ) {
			cendata_t *pt = &(data->at(j));
			if( pt->x > nct[i-1] ) {
				r += weights->at(j);
				if( pt->x <= nct[i] && pt->e == 1 ) d += weights->at(j);
			}
		}
		double p = (r==0)?1:(1-d/r);
		surv_prob[i] = surv_prob[i-1] * p;
	}
	surv_prob.push_back(0.0);
	nct.push_back(INFINITY);

	// set to member variables
	mProb = surv_prob;
	mTime = nct;

	// convert survival probability to density
	for( unsigned int i=0 ; i<nct.size()-1 ; i++ ) {
		region_t region_one = {nct[i],nct[i+1]};
		Region region;
		region.push_back(region_one);

		deninfo_t deninfo;
		deninfo.prob = surv_prob[i]-surv_prob[i+1];
		deninfo.den = deninfo.prob/(nct[i+1]-nct[i]);

		mWKMDensity.addDensity(&region,&deninfo);
	}

	return;
}

double WKM::getSurvProb(double t)
{
	if( t < mTime.front() ) return 1.0;
	if( isinf(t) ) return 0.0;
	if( t >= mTime.back() ) return mProb.back();
	unsigned int idx = 0;
	for( unsigned int i=0 ; i<mTime.size() ; i++ ) {
		if( t == mTime[i] ) { idx = i; break; } 
		else if( t < mTime[i] ) { idx = i-1; break; }
	}
	return mProb[idx];
}

double WKM::getCondSurvProb(double t,double t0)
{
	if( t <= t0 ) return 1.0;
	double p1 = getSurvProb(t);
	double p2 = getSurvProb(t0);
	if( p2 == 0.0 || p1 == 0.0 ) return 0.0;
	if( p1 >= p2 ) return 1.0;	
	return( ROUND(p1/p2) );
}


