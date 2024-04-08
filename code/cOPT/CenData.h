#ifndef _CENDATA_H_
#define _CENDATA_H_

#include <vector>

using namespace std;

// for one sample point
typedef struct _cendata_t {
	double x;
	int e;
} cendata_t;

// multiple sample points in multidimension
class CenData
{
public:
	CenData();
	~CenData();

	int getDimension() { return mDimension; }
	int getDataSize() { return ((int)(mData[0].size())); }
	vector<cendata_t>* getData(int idx) { return( &(mData[idx]) ); }
	cendata_t* getDataSingle(int nidx,int idx) { return( &((mData[idx])[nidx]) ); }

	bool readDataFromFile(char *filename);
	bool writeDataToFile(char *filename);

private:
	int mDimension;
	vector<cendata_t> mData[16];
};

#endif


