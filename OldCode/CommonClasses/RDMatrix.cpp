#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "RDMatrix.hpp"
#include "CVector.hpp"
#include "Specs.hpp"
using namespace std;

//Default constructor.  Quick but not useable.
RDMatrix::RDMatrix(){
	xBins = 1;
	yBins = 1;
	fBins = 1;
	nElements = 1;
	s = NULL;
    entry = new double[1];
	entry[0] = 0;
}

//Creates a blank array to store to elements of a real, diagonal matrix with size from specs.txt
RDMatrix::RDMatrix(Specs *s){
	loadSpecs(s);
	entry = new double[nElements];
	for (int i = 0; i < nElements; i++){
		entry[i] = 0;
	}
}


//Creates a blank array to store the elements of a real, diagonal matrix
RDMatrix::RDMatrix(Specs *s, int N){
	loadSpecs(s);
	nElements = N;
    entry = new double[nElements];
	for (int i = 0; i < nElements; i++){
		entry[i] = 0;
	}
}

//Loads in the diagonal of a real, diagonal matrix from a file with size from specs.txt
RDMatrix::RDMatrix(Specs *s, string filename){
	loadSpecs(s);
	fstream infile(filename.c_str(),fstream::in);
	entry = new double[nElements];
    for (int n = 0; n < nElements; n++){
        infile >> entry[n];
    }
	infile.close();
}


//Copy constructor.
RDMatrix::RDMatrix(const RDMatrix& copy){
	xBins = copy.xBins;
	yBins = copy.yBins;
	fBins = copy.fBins;
	nElements = copy.nElements;
	s = copy.s;
	entry = new double[nElements];
	for (int i = 0; i < nElements; i++){
		entry[i] = copy.entry[i];
    }
}

//Assignment operator
RDMatrix& RDMatrix::operator=(const RDMatrix& x){
	if (this == &x) return *this; //prevents self-assignment
	delete[] entry;
	xBins = x.xBins;
	yBins = x.yBins;
	fBins = x.fBins;
	nElements = x.nElements;
	s = x.s;
    entry = new double[nElements];
	for (int i = 0; i < nElements; i++){
        entry[i] = x.entry[i];
    }
	return *this;
}


//Destructor
RDMatrix::~RDMatrix(){
    delete[] entry;
}


//Loads relevant configuration information from file
void RDMatrix::loadSpecs(Specs *s){	
	this->s = s;
	xBins = s->xBins;
	yBins = s->yBins;
	fBins = s->fBins;
	nElements = xBins*yBins*fBins;
}


//Performs matrix multiplication with this matrix on the left of the given vector, overwriting the *
CVector RDMatrix::operator*(const CVector& x){
    CVector result = CVector(s);
	for (int i = 0; i < nElements; i++){
        result.real[i] = x.real[i] * entry[i];
		result.imag[i] = x.imag[i] * entry[i];
    }
	return result;
}

//This function outputs the entirety of the matrix
void RDMatrix::printAll(){
    cout << "This is a " << nElements << " by " << nElements << " real, diagonal matrix." << endl;
    cout << "It's entries are:" << endl;
    for (int n = 0; n < nElements; n++){
        cout << "(" << n << "," << n << "): " << entry[n] << endl;
    }
}
