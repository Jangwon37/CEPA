#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "CenData.h"
#include "Region.h"
#include "WKM.h"
#include "OPT.h"
#include "PHI.h"
#include "MyMath.h"

double pINF = 1.0;

using namespace std;

void help()
{
	cerr << "copt [options] data_file" << endl;
	cerr << "v0.2, Mar. 7th, 2012, by Junhee Seok, Stanford University" << endl;
	cerr << "options" << endl;
	cerr << "   -o filename    : output density file (default screen out)" << endl;
	cerr << "   -s b-e,b-e,... : region of interests (default data range)" << endl;
	cerr << "   -p filename    : using prior density in filename (default not using prior)" << endl;
	cerr << "   -d integer     : maximum depth of OPT (default 8)" << endl;
	cerr << "   -n float       : minimum num of effective points per region (default 1)" << endl;
	cerr << "   -i integer     : maximum number of iterations (default 5)" << endl;
	cerr << "   -e float       : criteria to stop iterations (default 0.001)" << endl;
	cerr << "   -r float       : rho (default 0.5)" << endl;
	cerr << "   -a float       : alpha (default 0.5)" << endl;
	cerr << "   -P filename    : output phi table (default none)" << endl;
	cerr << "   -O             : perform OPT once (default iteration)" << endl;
	cerr << "   -I             : partitioning at integer points (default FALSE)" << endl;
	cerr << "   -M             : mixed mode on (default off)" << endl;
	cerr << "   -v             : verb mode" << endl;
	cerr << "   -h             : show this message" << endl;
}

bool parse_roi(char *_str,Region *roi)
{
	char str[10240];
	strcpy(str,_str);
	char *p[32];
	int np = 1, n = strlen(str);
	p[0] = str;
	for( int i=0 ; i<n ; i++ ) {
		if( str[i] == '-' || str[i] == ',' ) { str[i] = '\0'; p[np] = &str[i+1]; np++; }
	}
	if( np%2 > 0 ) return false;
	int nd = np/2;

	region_t r;
	for( int i=0 ; i<nd ; i++ ) {
		r.begin = atof(p[2*i]);
		r.end = atof(p[2*i+1]);
		roi->push_back(r);
	}

	return true;
}

int main(int argc, char *argv[])
{
	// default parameters
	char *in_file = NULL;
	char *out_file = NULL;
	char *phi_file = NULL;
	char *prior_file = NULL;
	char *roi_str = NULL;
	int max_depth = 8;
	double min_num_points = 1.0;
	int max_iter = 5;
	double min_error = 0.001;
	double rho = 0.5;
	double alpha = 0.5;
	bool int_part = false;
	bool do_once = false;
	bool mixed_mode = false;
	bool verb_mode = false;

	// read parameters
	int c;
	opterr = 0;
	while( (c=getopt(argc,argv,"P:p:o:s:d:n:i:e:r:a:OIMvh")) != -1 ) {
		switch(c) {
			case 'p': prior_file = optarg; break;
			case 'o': out_file = optarg; break;
			case 's': roi_str = optarg; break;
			case 'd': max_depth = atoi(optarg); break;
			case 'n': min_num_points = atof(optarg); break;
			case 'i': max_iter = atoi(optarg); break;
			case 'e': min_error = atof(optarg); break;
			case 'r': rho = atof(optarg); break;
			case 'a': alpha = atof(optarg); break;
			case 'O': do_once = true; break;
			case 'I': int_part = true; break;
			case 'P': phi_file = optarg; break;
			case 'M': mixed_mode = true; break;
			case 'h': help(); exit(0); break;
			case 'v': verb_mode = true; break;
		}
	}
	if( optind == argc ) { help(); exit(0); }
	in_file = argv[optind];

	// read data
	CenData Data;
	if( !Data.readDataFromFile(in_file) ) { exit(0); }

	// region of interests
	Region ROI, *pROI = NULL;
	if( roi_str != NULL ) {
		if( !parse_roi(roi_str,&ROI) ) { 
			cerr << "Invalid ROI format" << endl;
			help(); exit(0); 
		}
		if( (int)(ROI.size()) != Data.getDimension() ) { 
			cerr << "The dimension of the given ROI does not agree with the data dimension" << endl;
			exit(0);
		}
		pROI = &ROI;
	}

	// OPT
	OPT OPT(&Data,pROI);
	OPT.setMaxDepth(max_depth);
	OPT.setMinNumPoints(min_num_points);
	OPT.setMaxIter(max_iter);
	OPT.setMinErr(min_error);
	OPT.setRho(rho);
	OPT.setAlpha(alpha);
	OPT.setIntPart(int_part);
	OPT.setMixedMode(mixed_mode);
	if( verb_mode ) OPT.printOptOptions();

	// prior density
	Density Sp, *pSp = NULL;
	if( prior_file && do_once ) {
		Sp.readDensityFromFile(prior_file);
		pSp = &Sp;
	}

	// run opt
	if( do_once ) OPT.runOptOne(pSp);
	else OPT.runOpt();

	// output
	Density *Den = OPT.getDensity();
	if( out_file ) Den->writeDensityToFile(out_file);
	else Den->printDensity();

	// optional PHI table output
	PHI *Phi = OPT.getPHI();
	if( phi_file ) Phi->writePhiToFile(phi_file);

	return 0;
}



