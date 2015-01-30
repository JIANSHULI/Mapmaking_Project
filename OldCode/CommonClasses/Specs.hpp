#include <fstream>

#ifndef SPECIFICATIONS
#define SPECIFICATIONS

using namespace std;

class Specs {

    public:
		Specs();
		Specs(string filename);
		~Specs();
		int xBins, yBins, fBins, zeroPad, FisherMCNum;
		double CGBound, xyLength, fLength, fStart, PreconEVThreshold;
		
};

#endif
