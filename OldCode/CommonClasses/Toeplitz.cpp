#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "Toeplitz.hpp"
#include "CVector.hpp"
#include "Specs.hpp"
#include "fftw3.h"
#include <string>
#include <math.h>

using namespace std;

const double pi = 3.14159265;


//Default constructor.  Quick but not useable.
Toeplitz::Toeplitz(){
	nElements = 1;
	xDim = false;
	yDim = false;
	fDim = false;
	xPos = 0;
	yPos = 0;
    entry = new double[1];
	entry[0] = 0;
}

//Loads in one dimension of the toeplitz matrices that makeup U and G
Toeplitz::Toeplitz(int n, string filename, string dim){
	xDim = false; yDim = false; fDim = false; xPos = -1; yPos = -1;
	nElements = n;
	if (dim.compare("x") == 0){
		xDim = true;
	} else if (dim.compare("y") == 0){
		yDim = true;
	} else if (dim.compare("f") == 0 || dim.compare("z") == 0){
		fDim = true;
	} else cout << "Warning: no toeplitz dimension specified!" << endl;
	fstream infile(filename.c_str(),fstream::in);
	entry = new double[nElements];
    for (int n = 0; n < nElements; n++){
        infile >> entry[n];
    }
	infile.close();
}

//Loads in a toeplitz matrix to represent one source in R
Toeplitz::Toeplitz(int n, string filename, int x, int y){
	xDim = false; yDim = false; fDim = false;
	nElements = n;
	xPos = x;
	yPos = y;
	fstream infile(filename.c_str(),fstream::in);
	entry = new double[nElements];
    for (int n = 0; n < nElements; n++){
        infile >> entry[n];
    }
	infile.close();
}


//Copy constructor.
Toeplitz::Toeplitz(const Toeplitz& copy){
	nElements = copy.nElements;
	xDim = copy.xDim;
	yDim = copy.yDim;
	fDim = copy.fDim;
	xPos = copy.xPos;
	yPos = copy.yPos;
    entry = new double[nElements];
	for (int i = 0; i < nElements; i++){
		entry[i] = copy.entry[i];
    }
}

//Assignment operator
Toeplitz& Toeplitz::operator=(const Toeplitz& x){
	if (this == &x) return *this; //prevents self-assignment

	delete[] entry;
	nElements = x.nElements;
	xDim = x.xDim;
	yDim = x.yDim;
	fDim = x.fDim;
	xPos = x.xPos;
	yPos = x.yPos;
    entry = new double[nElements];
	for (int i = 0; i < nElements; i++){
        entry[i] = x.entry[i];
    }
	return *this;
}

//Destructor
Toeplitz::~Toeplitz(){
    delete[] entry;
}

//Adds together two toeplitz matrices
Toeplitz Toeplitz::operator+(const Toeplitz& T){
	Toeplitz Sum = T;
	for (int i = 0; i < nElements; i++) Sum.entry[i] += entry[i];
	return Sum;
}


//Performs matrix multiplication with this matrix on the left of the given vector, overwriting the *
CVector Toeplitz::operator*(const CVector& x){
	CVector product = x;
	int aMin = 0, aMax, bMin = 0, bMax, subBins;
	if (xDim) {
		aMax = x.yBins;
		bMax = x.fBins;
		subBins = x.xBins;
	} else if (yDim) {
		aMax = x.xBins;
		bMax = x.fBins;
		subBins = x.yBins;
	} else if (fDim) {
		aMax = x.xBins;
		bMax = x.yBins;
		subBins = x.fBins;
	} else {
		aMin = xPos;
		aMax = xPos + 1;
		bMin = yPos;
		bMax = yPos + 1;
		subBins = x.fBins;
	}
	for (int a = aMin; a < aMax; a++){
		for (int b = bMin; b < bMax; b++){
			double subReal[subBins], subImag[subBins];
			//Pull out the particular subvector corresponding to the correct dimension
			for (int i = 0; i < subBins; i++){
				if (xDim) {
					subReal[i] = x.real[b + a*x.fBins + i*x.fBins*x.yBins];
					subImag[i] = x.imag[b + a*x.fBins + i*x.fBins*x.yBins];
				} else if (yDim) {
					subReal[i] = x.real[b + i*x.fBins + a*x.fBins*x.yBins];
					subImag[i] = x.imag[b + i*x.fBins + a*x.fBins*x.yBins];
				} else {
					subReal[i] = x.real[i + b*x.fBins + a*x.fBins*x.yBins];
					subImag[i] = x.imag[i + b*x.fBins + a*x.fBins*x.yBins];
				}
			}
			
			//embed the Teoplitz matrix in a circulent matrix
			double circulant[subBins*2];
			for (int i = 0; i < subBins; i++) circulant[i] = entry[i];
			circulant[subBins] = 0;
			for (int i = 0; i < subBins - 1; i++) circulant[i + subBins + 1] = entry[subBins - 1 - i];
			if (subBins != nElements) cout << "Warning: dimension size and Toeplitz size do not match!" << endl;
			
			//zero-pad the input subvector
			double yReal[subBins*2], yImag[subBins*2];
			for (int i = 0; i < subBins*2; i++){
				if (i < subBins){
					yReal[i] = subReal[i];
					yImag[i] = subImag[i];
				} else {
					yReal[i] = 0;
					yImag[i] = 0;
				}
			}
			
			//Fourier Transform both vectors
			fftw_complex *circIn, *circOut, *yIn, *yOut, *productIn, *productOut;
			circIn = (fftw_complex*) fftw_malloc(subBins*2 * sizeof(fftw_complex));
			circOut = (fftw_complex*) fftw_malloc(subBins*2 * sizeof(fftw_complex));
			yIn = (fftw_complex*) fftw_malloc(subBins*2 * sizeof(fftw_complex));
			yOut = (fftw_complex*) fftw_malloc(subBins*2 * sizeof(fftw_complex));
			productIn = (fftw_complex*) fftw_malloc(subBins*2 * sizeof(fftw_complex));
			productOut = (fftw_complex*) fftw_malloc(subBins*2 * sizeof(fftw_complex));
			
			for (int i = 0; i < subBins*2; i++){
				circIn[i][0] = circulant[i];
				circIn[i][1] = 0;
				yIn[i][0] = yReal[i];
				yIn[i][1] = yImag[i];
			}
			fftw_plan circFT = fftw_plan_dft_1d(subBins*2, circIn, circOut, FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_plan yFT = fftw_plan_dft_1d(subBins*2, yIn, yOut, FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_execute(circFT);
			fftw_execute(yFT);
			
			//element by element complex multiplication
			for (int i = 0; i < subBins*2; i++){
				if (fabs(circOut[i][1]) > 1e-10) cout << "Warning: unexpected imaginary eigenvalues of C!" << endl;
				productIn[i][0] = circOut[i][0] * yOut[i][0];
				productIn[i][1] = circOut[i][0] * yOut[i][1];
			}
			
			//Fourier transform back the product to get the convolution
			fftw_plan productIFT = fftw_plan_dft_1d(subBins*2, productIn, productOut, FFTW_BACKWARD, FFTW_ESTIMATE);
			fftw_execute(productIFT);
			
			//Put back into original vector			
			for (int i = 0; i < subBins; i++){
				if (xDim) {
					product.real[b + a*x.fBins + i*x.fBins*x.yBins] = productOut[i][0]/subBins/2;
					product.imag[b + a*x.fBins + i*x.fBins*x.yBins] = 0;
				} else if (yDim) {
					product.real[b + i*x.fBins + a*x.fBins*x.yBins] = productOut[i][0]/subBins/2;
					product.imag[b + i*x.fBins + a*x.fBins*x.yBins] = 0;
				} else {
					product.real[i + b*x.fBins + a*x.fBins*x.yBins] = productOut[i][0]/subBins/2;
					product.imag[i + b*x.fBins + a*x.fBins*x.yBins] = 0;
				}
			}
			
			//Deallocate memory;
			fftw_destroy_plan(circFT); fftw_destroy_plan(yFT); fftw_destroy_plan(productIFT);
			fftw_free(circIn); fftw_free(circOut); fftw_free(yIn); fftw_free(yOut); fftw_free(productIn); fftw_free(productOut);
		}
	}
	return product;
}

vector< vector<double> > Toeplitz::gaussRandField2D(Toeplitz Y){ //TODO: fix this based on what I learned when preparing grad lunch
	if (!xDim) cout << "Warning: Incorrect Dimension Applied to 2D Gaussian Random Field" << endl;
	vector< vector<double> > gRandField(nElements, vector<double>(Y.nElements,0));
	
	//Circulant embedding
	int Nx = nElements*2-1;
	int Ny = Y.nElements*2-1;
	double cx[Nx], cy[Ny];
	for (int n = 0; n < nElements; n++) cx[n] = entry[n];
	for (int n = 0; n < Y.nElements; n++) cy[n] = Y.entry[n];
	for (int n = 0; n < nElements - 1; n++) cx[n + nElements] = entry[nElements - 1 - n];
	for (int n = 0; n < Y.nElements - 1; n++) cy[n + Y.nElements] = Y.entry[Y.nElements - 1 - n];
	
	//Fourier transform, square root, and multiply by gaussian random numbers
	fftw_complex *circIn, *circOut, *productIn, *productOut, *gaussIn, *gaussOut;
	circIn = (fftw_complex*) fftw_malloc(Nx*Ny * sizeof(fftw_complex));
	circOut = (fftw_complex*) fftw_malloc(Nx*Ny * sizeof(fftw_complex));
	productIn = (fftw_complex*) fftw_malloc(Nx*Ny * sizeof(fftw_complex));
	productOut = (fftw_complex*) fftw_malloc(Nx*Ny * sizeof(fftw_complex));
	gaussIn = (fftw_complex*) fftw_malloc(Nx*Ny * sizeof(fftw_complex));
	gaussOut = (fftw_complex*) fftw_malloc(Nx*Ny * sizeof(fftw_complex));
	for (int i = 0; i < Nx; i++){
		for (int j = 0; j < Ny; j++){
			circIn[i*Ny + j][0] = cx[i]*cy[j];
			circIn[i*Ny + j][1] = 0;
			double r1 = 1.0*rand()/RAND_MAX;
			double r2 = 1.0*rand()/RAND_MAX;
			gaussIn[i*Ny+j][0] = sqrt(-2.0*log(r1))*cos(2.0*pi*r2);
			gaussIn[i*Ny+j][1] = 0;
		}
	}
	fftw_plan circ2DFFT = fftw_plan_dft_2d(Nx, Ny, circIn, circOut, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan gauss2DFFT = fftw_plan_dft_2d(Nx, Ny, gaussIn, gaussOut, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(circ2DFFT);
	fftw_execute(gauss2DFFT);
	
	for (int i = 0; i < Nx; i++){
		for (int j = 0; j < Ny; j++){
			productIn[i*Ny + j][0] = sqrt(fabs(circOut[i*Ny + j][0])) * gaussOut[i*Ny + j][0];
			productIn[i*Ny + j][1] = sqrt(fabs(circOut[i*Ny + j][0])) * gaussOut[i*Ny + j][1];
		}
	}
	
	
	//Fourier transform back and normalize
	fftw_plan product2DFFT = fftw_plan_dft_2d(Nx, Ny, productIn, productOut, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(product2DFFT);
	for (int i = 0; i < nElements; i++){
		for (int j = 0; j < Y.nElements; j++){
			gRandField[i][j] = productOut[i*Ny + j][0] / pow(Nx*Ny,.5) / sqrt(1.0*Nx*Ny);
		}
	}
	
	//Deallocate memory;
	fftw_destroy_plan(circ2DFFT); fftw_destroy_plan(gauss2DFFT); fftw_destroy_plan(product2DFFT);
	fftw_free(circIn); fftw_free(circOut); fftw_free(productIn); fftw_free(productOut); fftw_free(gaussIn); fftw_free(gaussOut);

	return gRandField;
}

//This function outputs the entirety of the matrix
void Toeplitz::printAll(){
    cout << "This is a " << nElements << " by " << nElements << " symmetric Toeplitz matrix." << endl;
	if (xDim) cout << "It is for the x dimension." << endl;
	else if (yDim) cout << "It is for the y dimension." << endl;
	else if (fDim) cout << "It is for the f dimension." << endl;
	else cout << "It is for position (" << xPos << ", " << yPos << ")." << endl;
    cout << "Its first column is:" << endl;
    for (int n = 0; n < nElements; n++){
        cout << "(" << n << "," << n << "): " << entry[n] << endl;
    }
}

