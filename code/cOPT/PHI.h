#ifndef _PHI_H_
#define _PHI_H_

#include <map>
#include "Region.h"

typedef struct _phiinfo_t {
	double lphi;
	double lphi0;
	double lphia;
	double lphias[16];
	double den;
	double count;
	int depth;

	Region region;
	
	struct _phiinfo_t *parent;
	struct _phiinfo_t *child[16][2];
	int hits;
} phiinfo_t;

class PHI
{
public:
	PHI() {}
	~PHI() {}

	void setDimension(int d) { mDimension=d; }
	bool addPhi(Region *region,phiinfo_t *phiinfo);
	void clearPhi() { mPHI.clear(); }

	phiinfo_t* findPhi(Region *region);
	map<Region,phiinfo_t,RegionCompare>* getPhi() { return &mPHI; }

	bool writePhiToFile(char *filename);	
	void printPhi();
	
private:
	int mDimension;
	map<Region,phiinfo_t,RegionCompare> mPHI;
};

#endif

