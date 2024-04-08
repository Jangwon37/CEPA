#ifndef _REGION_H_
#define _REGION_H_

#include <vector>
#include <math.h>
#include "MyMath.h"

using namespace std;

typedef struct _region_t {
	double begin;
	double end;
} region_t;

typedef vector<region_t> Region;

class RegionCompare 
{
public:
	bool operator()(const Region a,const Region b) const {
		for( unsigned int i=0 ; i<a.size() ; i++ ) {
			if( a[i].begin < b[i].begin ) return true;
			else if( a[i].begin > b[i].begin ) return false;
			else if( a[i].end < b[i].end ) return true;
			else if( a[i].end > b[i].end ) return false;
		}
		return false;
	}
};

#endif

