#ifndef _DENSITY_H_
#define _DENSITY_H_

#include <map>
#include "Region.h"

typedef struct _deninfo_t {
	double prob;
	double den;
} deninfo_t;

class Density 
{
public:
	Density(int dim);
	Density() { mDimension=1; }
	~Density() {}

	bool addDensity(Region *region,deninfo_t *deninfo);
	void clearDensity() { mDensity.clear(); }
	map<Region,deninfo_t,RegionCompare>* getDensity() { return &mDensity; }

	double getProb(Region *region);
	double getCondProb(Region *region,Region *given);
	vector<double> getProbGrid(Region *region,int ngrid);
	void renormalize();

	bool readDensityFromFile(char *filename);
	bool writeDensityToFile(char *filename);	
	void printDensity();
	
private:
	double getArea(const Region *R,bool rep_inf);
	int mDimension;
	map<Region,deninfo_t,RegionCompare> mDensity;
};

#endif

