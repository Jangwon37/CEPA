/***********************************************************************************************
 *
 * Density Class
 *
 ***********************************************************************************************/

/***********************************************************************************************
 * Including Headers
 ***********************************************************************************************/

#include <iostream>
#include <fstream>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "Density.h"
#include "Region.h"
#include "MyMath.h"

using namespace std;

/***********************************************************************************************
 * Basic Functions
 ***********************************************************************************************/

Density::Density(int den)
{
	mDimension = den;
}

bool Density::addDensity(Region *region,deninfo_t *deninfo)
{
	mDensity.insert(make_pair(*region,*deninfo));
	return true;
}

/***********************************************************************************************
 * Probability Handlers
 ***********************************************************************************************/

double Density::getProb(Region *region)
{
	double prob = 0.0;

	for( map<Region,deninfo_t,RegionCompare>::iterator it=mDensity.begin() ; it!=mDensity.end() ; it++ ) {
		const Region* denreg = &(it->first);
		deninfo_t *deninfo = &(it->second);
		bool skip = false;

		// calculate an overlapped region
		double region_area = 1.0, overlapped_area = 1.0;
		for( unsigned int j=0 ; j<denreg->size() ; j++ ) {
			double b1 = region->at(j).begin;
			double b2 = denreg->at(j).begin;
			double e1 = region->at(j).end;
			double e2 = denreg->at(j).end;

			double b = MAX(b1,b2);
			double e = MIN(e1,e2);
			if( b >= e ) { skip = true; break; }

			if( isinf(e) ) { e = e2 = pINF; }
			region_area *= e2-b2;
			overlapped_area *= e-b;
		}
		double prob_one = 0.0;
		if( !skip & !isinf(region_area) ) { prob_one = deninfo->prob / region_area * overlapped_area; }
		prob += prob_one;
	}

	return prob;
}

double Density::getCondProb(Region *region,Region *given)
{
	int nd = region->size();

	// get overlapped region
	Region ov = *region;
	for( int i=0 ; i<nd ; i++ ) {
		double b = MAX(region->at(i).begin,given->at(i).begin);
		double e = MIN(region->at(i).end,given->at(i).end);
		if( b >= e ) return 0.0;	// no overlapping
		ov[i].begin = b;
		ov[i].end = e;
	}

	double p1 = getProb(&ov);
	double p2 = getProb(given);
	if( p1 == 0.0 || p2 == 0.0 ) return 0.0;
	if( p1 > p2 ) return 1.0;
	return ROUND(p1/p2);
}

vector<double> Density::getProbGrid(Region *region,int ngrid)
{
	int nd = region->size();

	int idx[16];
	memset(idx,0,sizeof(int)*16);

	int f[16];
	f[0] = 1;
	for( int i=1 ; i<=nd ; i++ ) f[i] = f[i-1] * (ngrid+1);
	int len = f[nd];

	vector<double> prob(len);
	Region R = *region;
	for( int i=0 ; i<len ; i++ ) {
		for( int j=0 ; j<nd ; j++ ) { idx[nd-j-1] = (i%f[j+1])/f[j]; }
		for( int j=0 ; j<nd ; j++ ) {
			double d = (region->at(j).end-region->at(j).begin)/ngrid;
			R[j].begin = region->at(j).begin + idx[j]*d;
			if( idx[j] == ngrid ) R[j].end = INFINITY;
			else R[j].end   = region->at(j).begin + (idx[j]+1)*d;
		}
		prob[i] = getProb(&R);
	}

	return prob;
}

/***********************************************************************************************
 * Utilities
 ***********************************************************************************************/

void Density::renormalize() 
{
	double s = 0.0;
	for( map<Region,deninfo_t,RegionCompare>::iterator it=mDensity.begin() ; it!=mDensity.end() ; it++ ) {
		deninfo_t *deninfo = &(it->second);
		s += deninfo->prob;
	}
	for( map<Region,deninfo_t,RegionCompare>::iterator it=mDensity.begin() ; it!=mDensity.end() ; it++ ) {
		const Region* denreg = &(it->first);
		deninfo_t *deninfo = &(it->second);
		deninfo->prob /= s;
		deninfo->den = deninfo->prob/getArea(denreg,false);
	}
}

double Density::getArea(const Region *R,bool rep_inf)
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

/***********************************************************************************************
 * Input/Output Functions
 ***********************************************************************************************/

bool Density::readDensityFromFile(char *filename)
{
	ifstream infile(filename);
	if( !infile.is_open() ) {
		cerr << "cannot open file: " << filename << endl;
		return false;
	}

	string line;
	int nline = 0, dim = 0;
	while( !infile.eof() ) {
		getline( infile, line );
		if( line.length() == 0 ) continue;
		if( line[line.length()-1] == '\r' ) { line.erase(line.length()-1) = '\0'; }
		nline++;
		if( line[0] == '#' ) continue;

		// fast reading
		char str[10240];
		assert(line.length() < sizeof(str));
		strcpy(str,line.c_str());
		char *p[32];
		int np = 1, n=strlen(str);
		p[0] = str;	
		for( int i=0 ; i<n ; i++ ) {
			if( str[i] == '\t' ) { str[i] = '\0'; p[np] = &str[i+1]; np++; }
		}

		// set fields
		dim = (np-2)/2;
		Region R;
		for( int i=0 ; i<dim ; i++ ) {
			region_t r;
			r.begin = atof(p[2*i]);	
			r.end = atof(p[2*i+1]);	
			R.push_back(r);
		}
		deninfo_t deninfo;
		deninfo.den = atof(p[np-2]);
		deninfo.prob = atof(p[np-1]);
		
		// put
		mDensity.insert(make_pair(R,deninfo));	
	}

	return true;
}

bool Density::writeDensityToFile(char *filename)
{
        ofstream outfile(filename);
        if( !outfile.is_open() ) {
                cerr << "cannot write file: " << filename << endl;
                return false;
        }
	for( map<Region,deninfo_t,RegionCompare>::iterator it=mDensity.begin() ; it!=mDensity.end() ; it++ ) {
		const Region* region = &(it->first);
		deninfo_t *deninfo = &(it->second);
		for( unsigned int i=0 ; i<region->size(); i++ ) {
			outfile << region->at(i).begin << "\t" << region->at(i).end << "\t";
		}
		outfile << deninfo->den << "\t" << deninfo->prob << endl;
	}
        outfile.close();

	return true;
}

void Density::printDensity()
{
	for( map<Region,deninfo_t,RegionCompare>::iterator it=mDensity.begin() ; it!=mDensity.end() ; it++ ) {
		const Region* region = &(it->first);
		deninfo_t *deninfo = &(it->second);
		for( unsigned int i=0 ; i<region->size(); i++ ) {
			cout << region->at(i).begin << "\t" << region->at(i).end << "\t";
		}
		cout << deninfo->den << "\t" << deninfo->prob << endl;
	}
}



