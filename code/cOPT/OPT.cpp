/***********************************************************************************************
 *
 * Main OPT Class
 *
 ***********************************************************************************************/


/***********************************************************************************************
 * Including Headers
 ***********************************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "MyMath.h"
#include "OPT.h"
#include "PHI.h"
#include "Density.h"

#define INDEN(x)	for( int _k=0; _k<2*((x)-1) ; _k++ ) cout << " ";

/***********************************************************************************************
 * Constructor and Deconstructor
 ***********************************************************************************************/

OPT::OPT(CenData* X, Region* A)
{
	mData = X;
	mDataSize = X->getDataSize();
	mDimension = mData->getDimension();
	mPHI.setDimension(mDimension);

	for( int i=0 ; i<mDimension ; i++ ) {
		region_t region;
		if( A == NULL ) {
			vector<cendata_t> *data = mData->getData(i);
			double min_x = INFINITY, max_x = -INFINITY;
			for( unsigned int j=0 ; j<data->size() ; j++ ) {
				min_x = MIN(data->at(j).x,min_x);
				max_x = MAX(data->at(j).x,max_x);
			}
			region.begin = ROUND(min_x*0.95);
			region.end = ROUND(max_x*1.05);
		} else {
			region = A->at(i);
		}
		mROI.push_back(region);
		region.end = INFINITY;
		mRT.push_back(region);
	}

	mMaxIter = 5;
	mMinErr = 0.001;
	mNumGrid = 10;
	mMaxDepth = 8;
	mMinNumPoints = 1.0;
	mRho = 0.5;
	mAlpha = 0.5;
	mIntPart = false;
	mMixedMode = false;

	mLogBetaForAlpha = lbeta(mAlpha,mAlpha);
	mLogRho = log(mRho);
	mLog1mRho = log(1.0-mRho);
	mLog1pDimension = log(1.0/mDimension);
	mRndCount = 0;

	mPHIRoot = NULL;
}

OPT::~OPT()
{
	clearPhiInfo(mPHIRoot);
	mPHIRoot = NULL;
}

/***********************************************************************************************
 * Main OPT Functions
 ***********************************************************************************************/

void OPT::runOpt()
{
	vector<double> prob;
	vector<double> pre_prob( (mNumGrid+1)*(mNumGrid+1),1.0/((mNumGrid+1)*(mNumGrid+1)) );
	Density Sp;

	int niter = 0; 
	double err = 1E10;
	while( niter < mMaxIter && err > mMinErr ) {
		runOptOne( (niter==0)?NULL:&Sp );
		prob = mOptDensity.getProbGrid(&mROI,mNumGrid);

		double s = 0.0;
		for( unsigned int i=0 ; i<prob.size() ; i++  ) {
			double d = prob[i] - pre_prob[i];
			s += d*d;
		}
		err = sqrt(s/prob.size()); 

		Sp = mOptDensity;
		pre_prob = prob;
		niter++;
	}
}

void OPT::runOptOne(Density* Sp)
{
	mPHI.clearPhi();
	mOptDensity.clearDensity();
	clearPhiInfo(mPHIRoot);
	mPHIRoot = NULL;

	vector<double> P(mDataSize,1.0);
	runOptOneSub(&mRT,&P,Sp,1,NULL);
	runOptPartition();	
}

phiinfo_t* OPT::runOptOneSub(Region* A,vector<double>* P,Density *Sp,int depth,phiinfo_t *phi_parent)
{
	/************************************************
	 * Enter 
	 ************************************************/

#if 0
	INDEN(depth); cout << "ENTER: ";
	for( unsigned int i=0 ; i<A->size() ; i++ ) { 
		cout << "[" << A->at(i).begin << "," << A->at(i).end << ")"; 
		if( i < A->size()-1 ) cout << " x ";
	}
	cout << ", depth=" << depth << ", neff=" << sumProb(P) << endl;
#endif
	
	// to be returned
	phiinfo_t *phiinfo = new phiinfo_t;
	phiinfo->parent = phi_parent;
	phiinfo->region = (*A);
	phiinfo->hits = 1;
	for( int i=0 ; i<16 ; i++ ) { phiinfo->child[i][0] = phiinfo->child[i][1] = NULL; }
	if( phi_parent == NULL ) {
		clearPhiInfo(mPHIRoot);
		mPHIRoot = phiinfo;
	}

	/************************************************
	 * Try to find
	 ************************************************/

#if 0
	phiinfo_t *pre_phi = mPHI.findPhi(A);
	if( pre_phi ) {
		(*phiinfo) = (*pre_phi);

#if 0
	INDEN(depth); cout << "EXIT : ";
	for( unsigned int i=0 ; i<A->size() ; i++ ) { 
		cout << "[" << A->at(i).begin << "," << A->at(i).end << ")"; 
		if( i < A->size()-1 ) cout << " x ";
	}
	cout << ", depth=" << depth << ", lphi=" << phiinfo->lphi << ", FOUND FROM TABLE" << endl;
#endif

		return phiinfo;
	}
#endif

	/************************************************
	 * Prepare
	 ************************************************/

	// area of this region
	double area = getArea(A,true);
	if( area < 1E-10 ) { cerr << "too small area: " << area <<  endl; }

	// number of effective data points in this region
	double nall = sum_vector(P);
	vector<double> neff(mDimension,0.0);
	for( int i=0 ; i<mDimension ; i++ ) {
		vector<cendata_t>* data = mData->getData(i);
		for( unsigned int j=0 ; j<data->size() ; j++ ) {
			if( 	A->at(i).begin <= data->at(j).x && data->at(j).x < A->at(i).end &&
				P->at(j) > 0 && data->at(j).e == 1 ) neff[i] += P->at(j);
		}
	}

	// partition points
	vector<double> part_pt(mDimension);
	for( int i=0 ; i<mDimension ; i++ ) {
		part_pt[i] = ( A->at(i).begin + MIN(A->at(i).end,mROI[i].end) )/2;
		part_pt[i] = ROUND(part_pt[i]);
		if( mIntPart ) part_pt[i] = round(part_pt[i]);
	}

	// check stop criteria
	bool stop_partition = false;
	if( depth > mMaxDepth ) stop_partition = true;
	if( nall <= mMinNumPoints ) stop_partition = true;
	bool alltrue = true;
	for( int i=0 ; i<mDimension ; i++ ) { if( neff[i] > mMinNumPoints ) { alltrue = false; break; } }
	if( alltrue ) stop_partition = true;
	for( int i=0 ; i<mDimension ; i++ ) {
		double lower = ROUND(A->at(i).begin);
		double upper = ROUND(MIN(A->at(i).end,mROI[i].end));
		if( part_pt[i] <= lower || part_pt[i] >= upper ) { stop_partition = true; break; }
	}

	// if stop partitioning and on edges, partition one more
	if( stop_partition ) {
		for( int i=0 ; i<mDimension ; i++ ) {
			if( A->at(i).begin != mROI[i].end && isinf(A->at(i).end) ) { part_pt[i] = mROI[i].end; stop_partition = false; } 
			else { part_pt[i] = -INFINITY; }
		}
	}

	/************************************************
	 * Partitioning
	 ************************************************/

	if( stop_partition ) {
		phiinfo->lphi0 = log(1/area)*nall;
		phiinfo->lphi = phiinfo->lphi0;	
		phiinfo->lphia = -INFINITY;
		for( int i=0 ; i<mDimension ; i++ ) phiinfo->lphias[i] = -INFINITY;
	} else {
		// try to partition in one direction
		for( int i=0 ; i<mDimension ; i++ ) {
			phiinfo->lphias[i] = -INFINITY;
			if( A->at(i).begin < part_pt[i] && part_pt[i] < A->at(i).end ) {
				Region a1=(*A), a2=(*A);
				a1[i].end = a2[i].begin = part_pt[i];
				vector<double> p1(mDataSize,0.0);
				vector<double> p2(mDataSize,0.0);
				if( neff[i] > 0 ) {
					propagate(A,&a1,&a2,P,Sp,i,&p1,&p2);
					for( unsigned int j=0 ; j<p1.size() ; j++ ) { p1[j] *= P->at(j); p2[j] *= P->at(j); }
				} else {
					p2 = (*P);
				}

				phiinfo_t* phiinfo1 = runOptOneSub(&a1,&p1,Sp,depth+1,phiinfo);
				phiinfo_t* phiinfo2 = runOptOneSub(&a2,&p2,Sp,depth+1,phiinfo);
				phiinfo->child[i][0] = phiinfo1;
				phiinfo->child[i][1] = phiinfo2;

				double lphia = 0.0;
				lphia += lbeta(phiinfo1->count+mAlpha,phiinfo2->count+mAlpha) - mLogBetaForAlpha;
				lphia += phiinfo1->lphi + phiinfo2->lphi;
				phiinfo->lphias[i] = lphia;
			}
		}
		// summarize
		phiinfo->lphi0 = log(1/area)*nall;
		phiinfo->lphia = -INFINITY;
		for( int i=0 ; i<mDimension ; i++ ) {
			phiinfo->lphia = add_logged_values(phiinfo->lphia,mLog1pDimension+phiinfo->lphias[i]);
		}
		phiinfo->lphi = add_logged_values( mLogRho+phiinfo->lphi0, mLog1mRho+phiinfo->lphia );
	}
	phiinfo->den = nall/((double)mDataSize)/area;
	phiinfo->count = nall;
	phiinfo->depth = depth;

	/************************************************
	 * Update PHI table
	 ************************************************/

	phiinfo_t *pphi = mPHI.findPhi(A);
	if( pphi == NULL ) {
		mPHI.addPhi(A,phiinfo);
	} else if( mMixedMode ) {
		double w = (double)(pphi->hits), tmp;
		tmp = ( pphi->lphi*w + phiinfo->lphi )/(w+1.0); pphi->lphi = tmp;
		tmp = ( pphi->lphi0*w + phiinfo->lphi0 )/(w+1.0); pphi->lphi0 = tmp;
		tmp = ( pphi->lphia*w + phiinfo->lphia )/(w+1.0); pphi->lphia = tmp;
		for( int i=0 ; i<mDimension ; i++ ) {
			tmp = ( pphi->lphias[i]*w + phiinfo->lphias[i] )/(w+1.0); pphi->lphias[i] = tmp;
		}
		tmp = ( pphi->den*w + phiinfo->den )/(w+1.0); pphi->den = tmp;
		tmp = ( pphi->count*w + phiinfo->count )/(w+1.0); pphi->count = tmp;
		pphi->hits++;
	}

	/************************************************
	 * Exit
	 ************************************************/

#if 0
	INDEN(depth); cout << "EXIT : ";
	for( unsigned int i=0 ; i<A->size() ; i++ ) { 
		cout << "[" << A->at(i).begin << "," << A->at(i).end << ")"; 
		if( i < A->size()-1 ) cout << " x ";
	}
	cout << ", depth=" << depth << ", lphi=" << phiinfo->lphi << " :: ";
	cout << phiinfo->lphi0 << " " << phiinfo->lphia << " " << phiinfo->lphias[0] << " " << phiinfo->lphias[1] << endl;
#endif

	return( phiinfo );
}

/***********************************************************************************************
 * OPT Partitioning Functions
 ***********************************************************************************************/

void OPT::runOptPartition()
{
	runOptPartitionSub(mPHIRoot,1);
	mOptDensity.renormalize();
}

void OPT::runOptPartitionSub(phiinfo_t *phi,int depth)
{
	/************************************************
	 * Parameters
	 ************************************************/

	if( phi == NULL ) return;
	Region* A = &(phi->region);

	/************************************************
	 * Enter
	 ************************************************/

#if 0
	INDEN(depth); cout << "ENTER: ";
	for( unsigned int i=0 ; i<A->size() ; i++ ) { 
		cout << "[" << A->at(i).begin << "," << A->at(i).end << ")"; 
		if( i < A->size()-1 ) cout << " x ";
	}
	cout << ", depth=" << depth << endl; 
#endif

	/************************************************
	 * Finding PHI from table
	 ************************************************/

	if( mMixedMode ) {
		phiinfo_t *pphi = mPHI.findPhi(A);
		if( pphi ) {
			phi->lphi = pphi->lphi;
			phi->lphi0 = pphi->lphi0;
			phi->lphia = pphi->lphia;
			for( int i=0 ; i<mDimension ; i++ ) phi->lphias[i] = pphi->lphias[i];
			phi->den = pphi->den;
			phi->count = pphi->count;
		}
	}

	/************************************************
	 * Prepare
	 ************************************************/

	// density info to be added
	deninfo_t deninfo;
	deninfo.prob = phi->count/mDataSize;
	deninfo.den = deninfo.prob/getArea(A,false);

	bool stop_partition = false;
	bool less_lphia = ( (mLogRho+phi->lphi0) > (mLog1mRho+phi->lphia) );
	if( less_lphia ) stop_partition = true;
	bool all_child_null = true;
	for( int i=0 ; i<mDimension ; i++ ) {
		if( phi->child[i][0] != NULL ) { all_child_null = false; break; }
	}
	if( all_child_null ) stop_partition = true;

	if( stop_partition & !all_child_null ) {
		bool any_open_space = false;
		for( int i=0 ; i<mDimension ; i++ ) {
			if( isinf(phi->region[i].end) ) { any_open_space = true; break; }
		}
		if( any_open_space ) stop_partition = false;
	}

	/************************************************
	 * Partitioning
	 ************************************************/

	if( stop_partition ) {
		mOptDensity.addDensity(A,&deninfo);
	} else {
		double max_lphia = -INFINITY; 
		int max_lphia_idx = -1;
		for( int i=0 ; i<mDimension ; i++ ) {
			if( max_lphia < phi->lphias[i] ) {
				max_lphia = phi->lphias[i];
				max_lphia_idx = i;
			} else if( max_lphia == phi->lphias[i] ) {
				if( mRndCount++ % mDimension == i ) {
					max_lphia = phi->lphias[i];
					max_lphia_idx = i;
				}
			}
		}
		if( max_lphia_idx < 0 ) {
			mOptDensity.addDensity(A,&deninfo);
		} else {
			runOptPartitionSub(phi->child[max_lphia_idx][0],depth+1);
			runOptPartitionSub(phi->child[max_lphia_idx][1],depth+1);
		}
	}

	/************************************************
	 * Exit
	 ************************************************/

#if 0
	INDEN(depth); cout << "EXIT : ";
	for( unsigned int i=0 ; i<A->size() ; i++ ) { 
		cout << "[" << A->at(i).begin << "," << A->at(i).end << ")"; 
		if( i < A->size()-1 ) cout << " x ";
	}
	cout << ", depth=" << depth << ", prob=" << deninfo.prob << endl; 
#endif

	return;
}


/***********************************************************************************************
 * Probability Propagation
 ***********************************************************************************************/

// calculate the conditional probability in A1 given in A=A1+A2
void OPT::propagate(Region *A,Region* A1,Region* A2,vector<double>* Pin,Density* Sp,int idx,vector<double>* Pout1,vector<double>* Pout2)
{
	vector<cendata_t>* data = mData->getData(idx);
	region_t* a = &(A1->at(idx));

	// non-censored data
	for( int i=0 ; i<mDataSize ; i++ ) {
		if( a->begin <= data->at(i).x && data->at(i).x < a->end && data->at(i).e == 1 ) Pout1->at(i) = 1.0;
	}

	// censored data
	if( Sp == NULL ) {
		// given no prior distribution

		// weights
		vector<double> w(mDataSize,0.0);
		for( int i=0 ; i<mDataSize ; i++ ) {
			int setflag = 0;	// 0: no weight, 1: weight 1, 2: prob weight
			if( A->at(idx).begin <= data->at(i).x && data->at(i).x < A->at(idx).end ) {
				setflag = 1;
				for( int j=0 ; j<mDimension ; j++ ) {
					if( j == idx ) continue; // skip the target direction
					cendata_t* d = mData->getDataSingle(i,j);
					if( d->x >= A->at(j).end ) { setflag = 0; break; }	// no weight
					if( d->e == 0 || d->x <= A->at(j).begin ) { setflag = 2; }
				}
			}
			switch( setflag ) {
			case 1: w[i] = 1.0; break;
			case 2: w[i] = Pin->at(i); break;
			}
		}

		// collect data for local KM
		vector<cendata_t> local_data;
		vector<double> local_weight;
		for( int i=0 ; i<mDataSize ; i++ ) {
			if( w[i] > 0 ) {
				local_data.push_back(data->at(i));
				local_weight.push_back(w[i]);
			}
		}

		// local KM
		if( local_data.size() > 0 ) {
			WKM wkm(&local_data,&local_weight);
			for( int i=0 ; i<mDataSize ; i++ ) {
				bool alltrue = false;
				if( data->at(i).e == 0 && data->at(i).x < A1->at(idx).end ) {
					alltrue = true;
					for( int j=0 ; j<mDimension ; j++ ) {
						if( j == idx ) continue;
						cendata_t *d = mData->getDataSingle(i,j);
						if( d->x >= A->at(j).end ) { alltrue=false; break; }
					}
				}
				if( alltrue ) {
					double p1 = wkm.getCondSurvProb(A1->at(idx).begin,data->at(i).x);
					double p2 = wkm.getCondSurvProb(A2->at(idx).begin,data->at(i).x);
					if( isnan(p1) || isnan(p2) ) {
						cerr << "t>" << A1->at(idx).begin << "\tt>" << data->at(i).x << "\t" << p1 << endl;
						cerr << "t>" << A2->at(idx).begin << "\tt>" << data->at(i).x << "\t" << p2 << endl;
					}
					Pout1->at(i) = p1-p2;
				}
			}
		}

	} else {
		// given prior distribution

		Region r = (*A), g = (*A);
		for( int i=0 ; i<mDataSize ; i++ ) {

			if( data->at(i).e == 1 ) continue;	// non-censored were already considered
			if( Pin->at(i) == 0 ) continue;		// we don't have to consider if prior prob is 0
			if( data->at(i).x >= a->end ) continue;
			bool in_range = true;
			for( int j=0 ; j<mDimension ; j++ ) { 
				cendata_t* d = mData->getDataSingle(i,j);	
				if( d->x >= A1->at(j).end ) { in_range=false; break; } 
			}
			if( !in_range ) continue;

			for( int j=0 ; j<mDimension ; j++ ) {
				cendata_t* d = mData->getDataSingle(i,j);	
				if( d->e == 0 ) {
					r[j].begin = g[j].begin = MAX(d->x,A1->at(j).begin);
					r[j].end = A1->at(j).end; g[j].end = A2->at(j).end;
				} else {
					r[j].begin = g[j].begin = d->x - DELTA;
					r[j].end   = g[j].end   = d->x + DELTA;
				}
			}
			Pout1->at(i) = Sp->getCondProb(&r,&g);
		}

	}

	// calculate the conditional prob. for A2
	for( int i=0 ; i<mDataSize ; i++ ) { Pout2->at(i) = 1-Pout1->at(i); }
}

/***********************************************************************************************
 * Utility Functions
 ***********************************************************************************************/

double OPT::getArea(Region *R,bool rep_inf)
{
	double area = 1.0;
	for( unsigned int i=0 ; i<R->size() ; i++ ) {
		if( isinf(R->at(i).end) & rep_inf ) {
			area *= pINF - R->at(i).begin;
		} else {
			area *= R->at(i).end - R->at(i).begin;
		}
	}
	return area;
}

double OPT::sumProb(vector<double>* P)
{
	double s = 0.0;
	for( unsigned int i=0 ; i<P->size() ; i++ ) s += P->at(i);
	return s;
}

void OPT::clearPhiInfo(phiinfo_t *phiinfo)
{
	if( phiinfo == NULL ) return;
	for( int i=0 ; i<16 ; i++ ) {
		clearPhiInfo(phiinfo->child[i][0]);
		clearPhiInfo(phiinfo->child[i][1]);
	}
	phiinfo = NULL;
}

/***********************************************************************************************
 * Misc. Functions
 ***********************************************************************************************/

void OPT::printOptOptions()
{
	cout << "OPT Options" << endl;
	cout << "Number of data points:                 " << mDataSize << endl;
	cout << "Demension:                             " << mDimension << endl;
	cout << "Overall region:                        ";
	for( unsigned int i=0 ; i<mRT.size() ; i++ ) { 
		cout << "[" << mRT[i].begin << "," << mRT[i].end << ")"; 
		if( i < mRT.size()-1 ) cout << " x ";
		else cout << endl;
	}
	cout << "Region of interest:                    ";
	for( unsigned int i=0 ; i<mROI.size() ; i++ ) { 
		cout << "[" << mROI[i].begin << "," << mROI[i].end << ")"; 
		if( i < mROI.size()-1 ) cout << " x ";
		else cout << endl;
	}
	cout << "Maximum iterations:                    " << mMaxIter << endl;
	cout << "Minimum error in iteration:            " << mMinErr << endl;
	cout << "Number of grids in error calculation:  " << mNumGrid << endl;
	cout << "Maximum depth:                         " << mMaxDepth << endl;
	cout << "Minimum number of points in a region:  " << mMinNumPoints << endl;
	cout << "Partition at integer points:           " << (mIntPart?"Yes":"No") << endl;
	cout << "rho:                                   " << mRho << endl;
	cout << "alpha:                                 " << mAlpha << endl;
}

