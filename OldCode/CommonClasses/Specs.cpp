#include "Specs.hpp"
using namespace std;

Specs::Specs(){
}

Specs::Specs(string filename){
	fstream infile(filename.c_str(),fstream::in);
	string dummy;
	infile >> dummy >> xBins;
	infile >> dummy >> yBins;
	infile >> dummy >> fBins;
	infile >> dummy >> zeroPad;
	infile >> dummy >> PreconEVThreshold;
	infile >> dummy >> CGBound;
	infile >> dummy >> FisherMCNum;
	infile >> dummy >> xyLength;
	infile >> dummy >> fLength;
	infile >> dummy >> fStart;
	infile.close();
}

Specs::~Specs(){
}
