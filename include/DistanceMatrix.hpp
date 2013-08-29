#ifndef CLASS_DISTANCE_MATRIX
#define CLASS_DISTANCE_MATRIX

#include <iostream>
#include <fstream>
//#include <stdlib.h>
//#include <vector>
#include <math.h>
//#include <stdexcept>

//#include <sstream> //debug only

#include "ProteinGroup.hpp"
#include "Alignment.hpp"



using namespace std;

class DistanceMatrix : public ProteinGroup {

private:
	
	
	enum {ala, arg, asn, asp, cys, gln, glu, gly, his, ileu, leu, lys, met, phe, pro, ser, thr, trp, tyr, val, SKIP_AA=100};
	

	static const double DISTANCE_NO_OVERLAP	= -1.0;
	static const double DISTANCE_INFINITY		= -1.0;
	static const double DISTANCE_TOO_SMALL		= .00001;

	vector<vector <double > > distances2D_;
		
	void init(){ProteinGroup::init(); distances2D_.clear();}
	static double compute_pmb_distance(vector <unsigned char> const& sequence1, vector <unsigned char> const& sequence2);
	static double compute_pmb_distance_quick(vector <unsigned char> const& sequence1, vector <unsigned char> const& sequence2);
	
	void read_from_filehandle(ifstream& inFile);
	void print_to_filehandle(ofstream& outFile) const;

public:
	//DistanceMatrix(){init();}
	//DistanceMatrix(string const& inFilepath) {init(); read_from_file(inFilepath);}
	//DistanceMatrix(Alignment alignment, bool const& quick) {init(); create_from_alignment(alignment, quick);}
	
	//mutators
	void read_from_file(string const& inFilepath);
	void create_from_alignment(Alignment alignment, bool const& quick);
	
	//non-mutators
	void print_to_file(string const& outFilepath) const;
	void get_distances_as_1d(vector <double>&) const;
};
//
#endif
