#ifndef _OPT_H_
#define _OPT_H_

#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#include "Region.h"
#include "WKM.h"
#include "Density.h"
#include "PHI.h"
#include "MyMath.h"

class OPT
{
public:
	OPT(CenData* X, Region* A);
	~OPT();

	// main functions
	void runOpt(); 	// iterative opt
	void runOptOne(Density* Sp);	// run opt once
	phiinfo_t* runOptOneSub(Region* A,vector<double>* P,Density *Sp,int depth,phiinfo_t *parent);	// recursive handler
	void runOptPartition();
	void runOptPartitionSub(phiinfo_t* phi,int depth);

	// output
	PHI* getPHI() { return &mPHI; }
	Density* getDensity() { return &mOptDensity; }

	// parameter settlers
	void setMaxIter(int p) { mMaxIter=p; }
	void setMinErr(double p) { mMinErr=p; }
	void setNumGrid(int p) { mNumGrid=p; }
	void setMaxDepth(int p) { mMaxDepth=p; }
	void setMinNumPoints(double p) { mMinNumPoints=p; }
	void setRho(double p) { mRho=p; mLogRho=log(mRho); mLog1mRho=log(1-mRho); }
	void setAlpha(double p) { mAlpha=p; mLogBetaForAlpha = lbeta(mAlpha,mAlpha); }
	void setIntPart(bool p) { mIntPart=p; }
	void setMixedMode(bool p) { mMixedMode=p; }

	// print
	void printOptOptions();

private:

	// parameters
	int mMaxIter;
	double mMinErr;
	int mNumGrid;
	int mMaxDepth;
	double mMinNumPoints;
	double mRho;
	double mAlpha;
	double mIntPart;
	double mMixedMode;

	// data
	int mDimension;
	Region mRT;	// total region
	Region mROI;	// region of interest
	CenData* mData;	// data
	int mDataSize;	// number of data points
	
	// output
	Density mOptDensity;
	PHI mPHI;
	phiinfo_t *mPHIRoot;

	// useful variables
	double mLogBetaForAlpha;
	double mLogRho, mLog1mRho;
	double mLog1pDimension;
	int mRndCount;

	// some functions
	void clearPhiInfo(phiinfo_t *phiinfo);
	double getArea(Region *R,bool rep_inf);
	double sumProb(vector<double>* P);
	void propagate(Region *A,Region* A1,Region* A2,vector<double>* Pin,Density* Sp,int idx,vector<double>* Pout1,vector<double>* Pout2);
};


#endif

