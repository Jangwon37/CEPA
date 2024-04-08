#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "CenData.h"
#include "MyMath.h"

using namespace std;

CenData::CenData() 
{
	mDimension = -1;
}

CenData::~CenData() 
{
}

bool CenData::readDataFromFile(char *filename)
{
	ifstream infile(filename);
	if( !infile.is_open() ) {
		cerr << "cannot open file: " << filename << endl;
		return false;
	}

	string line;
	int nline = 0;
	double max_value = -1.0;
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
		int ntokens = np;

		// check the number of items
		if( mDimension < 0 ) { 
			mDimension = ntokens/2; 
			if( ntokens%2 == 1 ) mDimension--;
			if( mDimension > 16 ) mDimension = 16;
		}
		int d = ntokens/2;
		if( ntokens%2 == 1 ) d--;
		if( mDimension == 0 || d != mDimension ) {
			cerr << "input file error: not enough fileds" << endl;
			cerr << filename << "(" << nline << "): " << line << endl;
			infile.close();
			return false;
		}

		// fill data
		for( int i=0 ; i<mDimension ; i++ ) {
			cendata_t cendata;
			cendata.x = ROUND(atof(p[2*i]));
			cendata.e = atoi(p[2*i+1]);
			if( cendata.e != 0 ) cendata.e = 1;
			mData[i].push_back(cendata);
			if( cendata.x > max_value ) max_value = cendata.x;
		}
	}

	infile.close();

	// set psuedo infinity
	pINF = round( 10.0 * max_value );

	// cout << "read " << nline << " samples from the file in " << mDimension << " dimension." << endl;

	return true;
}

bool CenData::writeDataToFile(char *filename) 
{
	ofstream outfile(filename);
	if( !outfile.is_open() ) {
		cerr << "cannot write file: " << filename << endl;
		return false;
	}

	for( unsigned int i=0 ; i<mData[0].size() ; i++ ) {
		for( int j=0 ; j<mDimension ; j++ ) {
			outfile << mData[j][i].x;
			if( mData[j][i].e == 0 ) outfile << "+";
			if( j != mDimension-1 )	outfile << "\t";
		}
		outfile << endl;
	}

	outfile.close();

	return true;
}

