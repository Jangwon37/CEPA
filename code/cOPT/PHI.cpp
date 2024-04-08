#include <iostream>
#include <fstream>
#include <string>
#include "PHI.h"

using namespace std;

bool PHI::addPhi(Region *region,phiinfo_t *phiinfo)
{
	mPHI.insert(make_pair(*region,*phiinfo));
	return true;
}

phiinfo_t* PHI::findPhi(Region *region)
{
	map<Region,phiinfo_t,RegionCompare>::iterator f = mPHI.find(*region);
	if( f == mPHI.end() ) return NULL;
	return &(f->second);
}

bool PHI::writePhiToFile(char *filename)
{
	ofstream outfile(filename);
	if( !outfile.is_open() ) {
		cerr << "cannot write file: " << filename << endl;
		return false;
	}
	for( map<Region,phiinfo_t,RegionCompare>::iterator it=mPHI.begin() ; it!=mPHI.end() ; it++ ) {
                const Region* region = &(it->first);
                phiinfo_t *phiinfo = &(it->second);
                for( unsigned int i=0 ; i<region->size(); i++ ) {
                        outfile << region->at(i).begin << "\t" << region->at(i).end << "\t";
                }
                outfile << phiinfo->lphi << "\t" << phiinfo->lphi0 << "\t" << phiinfo->lphia << "\t";
		for( int i=0 ; i<mDimension ; i++ ) outfile << phiinfo->lphias[i] << "\t";
		outfile << phiinfo->den << "\t" << phiinfo->count << "\t" << phiinfo->depth << endl;
        }
	outfile.close();

	return true;
}

void PHI::printPhi()
{
	for( map<Region,phiinfo_t,RegionCompare>::iterator it=mPHI.begin() ; it!=mPHI.end() ; it++ ) {
                const Region* region = &(it->first);
                phiinfo_t *phiinfo = &(it->second);
                for( unsigned int i=0 ; i<region->size(); i++ ) {
                        cout << region->at(i).begin << "\t" << region->at(i).end << "\t";
                }
                cout << phiinfo->lphi << "\t" << phiinfo->lphi0 << "\t" << phiinfo->lphia << "\t";
		for( int i=0 ; i<mDimension ; i++ ) cout << phiinfo->lphias[i] << "\t";
		cout << phiinfo->den << "\t" << phiinfo->count << "\t" << phiinfo->depth << endl;
        }

}

