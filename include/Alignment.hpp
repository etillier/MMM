#ifndef CLASS_ALIGNMENT
#define CLASS_ALIGNMENT

#include <cstring>
//#include <stdio.h>
#include <iostream>
#include <fstream>
//#include <stdlib.h>
//#include <string>
//#include <vector>
#include <map>
//#include <sstream>

#include "ProteinGroup.hpp"

#define BUFFER_SIZE 8192

using namespace std;

class Alignment : public ProteinGroup {

private:
	
	int nOfColumns_;
	int joinedAt_; //if the alignment is made by joining 2 alignments then it refers to the length of the first alignment
	
	vector <string> sequences_;
	
	void init(void){ProteinGroup::init(); nOfColumns_=-1; joinedAt_=0;}
	static int prepare_sequence(char *buffer);

public:
	Alignment(string const& inFilename, double gapThreshold) {init();read_alignment(inFilename, gapThreshold);}
	Alignment(string const& inFilename) {init();read_alignment(inFilename);}
	Alignment(int const joinedAt){init();joinedAt_ = joinedAt;}
	Alignment(){init();}
	~Alignment(){}
	
	void read_alignment(string const& inFilename, double gapThreshold);
	void read_alignment(string const& inFilename) {read_alignment(inFilename, -1);}
	
	void add_sequence(string header, string sequence);
	Alignment append_alignment(Alignment const& other);
		
	Alignment get_slices(vector<pair<int, int > > const& ranges) const;
	Alignment get_slice(int start, int end) const {pair <int,int> range (start,end); vector<pair<int, int > > ranges; ranges.push_back(range); get_slices(ranges);}
	void print_to_file(string const& outFilepath) const;

	int get_nOfColumns() const {validate(); return nOfColumns_;}
	int get_joinedAt() const {validate(); return joinedAt_;}
	string get_sequence(int index) const {validate(); return sequences_.at(index);}
};
//
#endif
