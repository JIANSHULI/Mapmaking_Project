#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "CVector.hpp"
#include "Specs.hpp"
#include "fftw3.h"

using namespace std;


//Default constructor.  Quick but not useable.
CVector::CVector(){
	nElements = 1;
	xBins = 1;
	yBins = 1;
	fBins = 1;
	zeroPad = 1;
    real = new double[1];
    imag = new double[1];
	s = NULL;
	real[0] = 0;
	imag[0] = 0;
}

//Creates two arrays initialized to all 0's of size determined by specs.txt (defaults to zeroPadded)
CVector::CVector(Specs *s){
    loadSpecs(s);
    real = new double[nElements];
    imag = new double[nElements];
    for (int i = 0; i < nElements; i++){
        real[i] = 0;
        imag[i] = 0;
    }
}


//Creates a complex vector with a specified number of elements
CVector::CVector(Specs *s, int N){
	loadSpecs(s);
	nElements = N;
    real = new double[nElements];
    imag = new double[nElements];
    for (int i = 0; i < nElements; i++){
        real[i] = 0;
        imag[i] = 0;
    }
}

//Creates a complex vector, but this time from data in a specially formatted file
CVector::CVector(Specs *s, string filename){
    loadSpecs(s);
	fstream infile(filename.c_str(),fstream::in);
	real = new double[nElements];
    imag = new double[nElements];
    for(int n = 0; n < nElements; n++){
        infile >> real[n] >> imag[n];
    }
	infile.close();
}

//Copy constructor
CVector::CVector(const CVector& copy){
	nElements = copy.nElements;
	xBins = copy.xBins;
	yBins = copy.yBins;
	fBins = copy.fBins;
	zeroPad = copy.zeroPad;
    real = new double[nElements];
    imag = new double[nElements];
	s = copy.s;
	for (int i = 0; i < nElements; i++){
        real[i] = copy.real[i];
        imag[i] = copy.imag[i];
    }
}

//Destructor
CVector::~CVector(){
    delete[] real;
    delete[] imag;
}

//Loads relevant configuration information from file
void CVector::loadSpecs(Specs *s){
	this->s = s;
	xBins = s->xBins;
	yBins = s->yBins;
	fBins = s->fBins;
	zeroPad = s->zeroPad;
	nElements = xBins*yBins*fBins;
}

//Assignment operator
CVector& CVector::operator=(const CVector& x){
	if (this == &x) return *this; //prevents self-assignment
	delete[] real;
	delete[] imag;
	s = x.s;
	loadSpecs(s);
	nElements = x.nElements;
	xBins = x.xBins;
	yBins = x.yBins;
	fBins = x.fBins;
	real = new double[nElements];
    imag = new double[nElements];
	for (int i = 0; i < nElements; i++){
        real[i] = x.real[i];
        imag[i] = x.imag[i];
    }
	return *this;
}

//returns the dot product of this object and the passed CVector as a length 1 CVector
CVector CVector::dot(const CVector& x){
	CVector dotProduct = CVector(s,1);
	for (int n = 0; n < nElements; n++){
        dotProduct.real[0] += x.real[n] * real[n] + x.imag[n] * imag[n];
        dotProduct.imag[0] += -x.real[n] * imag[n] + x.imag[n] * real[n];
    }
	return dotProduct;
}

//returns this object with every element multiplied by a complex scalar passed as a length 1 CVector;
CVector CVector::operator*(const CVector& c){
	CVector scalarMultiple = CVector(s);
	for (int n = 0; n < nElements; n++){
		scalarMultiple.real[n] = c.real[0] * real[n] - c.imag[0] * imag[n];
		scalarMultiple.imag[n] = c.real[0] * imag[n] + c.imag[0] * real[n];
	}
	return scalarMultiple;
}

//returns the sum of this object and the CVector passed
CVector CVector::operator+(const CVector& x){
	CVector sum = CVector(s);	
	for (int n = 0; n < nElements; n++){
		sum.real[n] = x.real[n] + real[n];
		sum.imag[n] = x.imag[n] + imag[n];
	}
	return sum;
}

//returns the difference of this object and the CVector passed
CVector CVector::operator-(const CVector& x){
	CVector difference = CVector(s);
	for (int n = 0; n < nElements; n++){
		difference.real[n] = real[n] - x.real[n];
		difference.imag[n] = imag[n] - x.imag[n];
	}
	return difference;
}

//returns the quotient of two complex numbers as a complex number (all represented as length 1 CVectors)
CVector CVector::operator/(const CVector& d){
	double divisor = d.real[0]*d.real[0] + d.imag[0]*d.imag[0];
	CVector quotient = CVector(s,1);
	quotient.real[0] = (d.real[0] * real[0] - d.imag[0] * imag[0])/divisor;
	quotient.imag[0] = (d.real[0] * imag[0] + d.imag[0] * real[0])/divisor;
	return quotient;
}

//returns true if the magnitude of each entry is less than a given number
bool CVector::allMagnitudesLessThan(double c){
	for (int n = 0; n < nElements; n++){
		if (sqrt(real[n]*real[n] + imag[n]*imag[n]) > c) return false;
	}
	return true;

}

//returns the sum of the magnitudes of all the entries
double CVector::magnitude(){
		double sum = 0;
		for (int j = 0; j < nElements; j++){
			sum += real[j]*real[j] + imag[j]*imag[j];
		}
		return sqrt(sum);
}

//Prints the contents of the vector to the console.
void CVector::printAll(string name){
    cout << endl << name << " is a vector with " << nElements << " element(s):" << endl;
    for (int i = 0; i < nElements; i++){
        cout << real[i] << " + " << imag[i] << "i" << endl;
    }
}

void CVector::printRealToFile(string filename){
	ofstream outfile;
	outfile.precision(20);	
	outfile.open(filename.c_str(), ios::trunc);	
	for (int i = 0; i < nElements; i++) outfile << real[i] <<  endl;
	outfile.close();
}

CVector CVector::ijk2uvw(){
	//Shift to start with the origin.
	CVector ijkShift = CVector(s);
	for (int n = 0; n < xBins*yBins; n++){
		for (int m = 0; m < fBins/2; m++){
			ijkShift.real[n*fBins + m] = real[n*fBins + fBins/2 + m];
			ijkShift.imag[n*fBins + m] = imag[n*fBins + fBins/2 + m];
			ijkShift.real[n*fBins + fBins/2 + m] = real[n*fBins + m];
			ijkShift.imag[n*fBins + fBins/2 + m] = imag[n*fBins + m];
		}
	}
	CVector ijkShiftTemp = CVector(s);
	for (int n = 0; n < xBins; n++){
		for (int m = 0; m < fBins*yBins/2; m++){
			ijkShiftTemp.real[n*fBins*yBins + m] = ijkShift.real[n*fBins*yBins + fBins*yBins/2 + m];
			ijkShiftTemp.imag[n*fBins*yBins + m] = ijkShift.imag[n*fBins*yBins + fBins*yBins/2 + m];
			ijkShiftTemp.real[n*fBins*yBins + fBins*yBins/2 + m] = ijkShift.real[n*fBins*yBins + m];
			ijkShiftTemp.imag[n*fBins*yBins + fBins*yBins/2 + m] = ijkShift.imag[n*fBins*yBins + m];
		}
	}
	for (int m = 0; m < fBins*yBins*xBins/2; m++){
		ijkShift.real[m] = ijkShiftTemp.real[fBins*yBins*xBins/2 + m];
		ijkShift.imag[m] = ijkShiftTemp.imag[fBins*yBins*xBins/2 + m];
		ijkShift.real[fBins*yBins*xBins/2 + m] = ijkShiftTemp.real[m];
		ijkShift.imag[fBins*yBins*xBins/2 + m] = ijkShiftTemp.imag[m];
	}
	
	//Perform the FFT
	fftw_complex *in, *out;
	in = (fftw_complex*) fftw_malloc(xBins*yBins*fBins * sizeof(fftw_complex));
	out = (fftw_complex*) fftw_malloc(xBins*yBins*fBins * sizeof(fftw_complex));
	for (int n = 0; n < nElements; n++){
		in[n][0] = ijkShift.real[n];
		in[n][1] = ijkShift.imag[n];
	}
	fftw_plan FT = fftw_plan_dft_3d(xBins, yBins, fBins, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(FT);
	CVector uvwShift = CVector(s);
	for (int n = 0; n < nElements; n++){
		uvwShift.real[n] = out[n][0];
		uvwShift.imag[n] = out[n][1];
	}
	fftw_destroy_plan(FT);
	fftw_free(in); fftw_free(out);
	
	//Shift back to zero frequency centered.
	CVector uvw = CVector(s);
	for (int m = 0; m < fBins*yBins*xBins/2; m++){
		uvw.real[m] = uvwShift.real[fBins*yBins*xBins/2 + m];
		uvw.imag[m] = uvwShift.imag[fBins*yBins*xBins/2 + m];
		uvw.real[fBins*yBins*xBins/2 + m] = uvwShift.real[m];
		uvw.imag[fBins*yBins*xBins/2 + m] = uvwShift.imag[m];
	}
	CVector uvwShiftTemp = CVector(s);
	for (int n = 0; n < xBins; n++){
		for (int m = 0; m < fBins*yBins/2; m++){
			uvwShiftTemp.real[n*fBins*yBins + m] = uvw.real[n*fBins*yBins + fBins*yBins/2 + m];
			uvwShiftTemp.imag[n*fBins*yBins + m] = uvw.imag[n*fBins*yBins + fBins*yBins/2 + m];
			uvwShiftTemp.real[n*fBins*yBins + fBins*yBins/2 + m] = uvw.real[n*fBins*yBins + m];
			uvwShiftTemp.imag[n*fBins*yBins + fBins*yBins/2 + m] = uvw.imag[n*fBins*yBins + m];
		}
	}
	for (int n = 0; n < xBins*yBins; n++){
		for (int m = 0; m < fBins/2; m++){
			uvw.real[n*fBins + m] = uvwShiftTemp.real[n*fBins + fBins/2 + m];
			uvw.imag[n*fBins + m] = uvwShiftTemp.imag[n*fBins + fBins/2 + m];
			uvw.real[n*fBins + fBins/2 + m] = uvwShiftTemp.real[n*fBins + m];
			uvw.imag[n*fBins + fBins/2 + m] = uvwShiftTemp.imag[n*fBins + m];
		}
	}
		
	return uvw;
}

CVector CVector::uvw2ijk(){

	//Shift to start with the origin.
	CVector uvwShift = CVector(s);
	for (int n = 0; n < xBins*yBins; n++){
		for (int m = 0; m < fBins/2; m++){
			uvwShift.real[n*fBins + m] = real[n*fBins + fBins/2 + m];
			uvwShift.imag[n*fBins + m] = imag[n*fBins + fBins/2 + m];
			uvwShift.real[n*fBins + fBins/2 + m] = real[n*fBins + m];
			uvwShift.imag[n*fBins + fBins/2 + m] = imag[n*fBins + m];
		}
	}
	CVector uvwShiftTemp = CVector(s);
	for (int n = 0; n < xBins; n++){
		for (int m = 0; m < fBins*yBins/2; m++){
			uvwShiftTemp.real[n*fBins*yBins + m] = uvwShift.real[n*fBins*yBins + fBins*yBins/2 + m];
			uvwShiftTemp.imag[n*fBins*yBins + m] = uvwShift.imag[n*fBins*yBins + fBins*yBins/2 + m];
			uvwShiftTemp.real[n*fBins*yBins + fBins*yBins/2 + m] = uvwShift.real[n*fBins*yBins + m];
			uvwShiftTemp.imag[n*fBins*yBins + fBins*yBins/2 + m] = uvwShift.imag[n*fBins*yBins + m];
		}
	}
	for (int m = 0; m < fBins*yBins*xBins/2; m++){
		uvwShift.real[m] = uvwShiftTemp.real[fBins*yBins*xBins/2 + m];
		uvwShift.imag[m] = uvwShiftTemp.imag[fBins*yBins*xBins/2 + m];
		uvwShift.real[fBins*yBins*xBins/2 + m] = uvwShiftTemp.real[m];
		uvwShift.imag[fBins*yBins*xBins/2 + m] = uvwShiftTemp.imag[m];
	}

	//Perform the FFT
	fftw_complex *in, *out;
	in = (fftw_complex*) fftw_malloc(xBins*yBins*fBins * sizeof(fftw_complex));
	out = (fftw_complex*) fftw_malloc(xBins*yBins*fBins * sizeof(fftw_complex));
	for (int n = 0; n < nElements; n++){
		in[n][0] = uvwShift.real[n];
		in[n][1] = uvwShift.imag[n];
	}
	fftw_plan FT = fftw_plan_dft_3d(xBins, yBins, fBins, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(FT);
	CVector ijkShift = CVector(s);
	for (int n = 0; n < nElements; n++){
		ijkShift.real[n] = out[n][0];
		ijkShift.imag[n] = out[n][1];
	}
	fftw_destroy_plan(FT);
	fftw_free(in); fftw_free(out);
	
	//Shift back to zero frequency centered.
	CVector ijk = CVector(s);
	for (int m = 0; m < fBins*yBins*xBins/2; m++){
		ijk.real[m] = ijkShift.real[fBins*yBins*xBins/2 + m];
		ijk.imag[m] = ijkShift.imag[fBins*yBins*xBins/2 + m];
		ijk.real[fBins*yBins*xBins/2 + m] = ijkShift.real[m];
		ijk.imag[fBins*yBins*xBins/2 + m] = ijkShift.imag[m];
	}
	CVector ijkShiftTemp = CVector(s);
	for (int n = 0; n < xBins; n++){
		for (int m = 0; m < fBins*yBins/2; m++){
			ijkShiftTemp.real[n*fBins*yBins + m] = ijk.real[n*fBins*yBins + fBins*yBins/2 + m];
			ijkShiftTemp.imag[n*fBins*yBins + m] = ijk.imag[n*fBins*yBins + fBins*yBins/2 + m];
			ijkShiftTemp.real[n*fBins*yBins + fBins*yBins/2 + m] = ijk.real[n*fBins*yBins + m];
			ijkShiftTemp.imag[n*fBins*yBins + fBins*yBins/2 + m] = ijk.imag[n*fBins*yBins + m];
		}
	}
	for (int n = 0; n < xBins*yBins; n++){
		for (int m = 0; m < fBins/2; m++){
			ijk.real[n*fBins + m] = ijkShiftTemp.real[n*fBins + fBins/2 + m];
			ijk.imag[n*fBins + m] = ijkShiftTemp.imag[n*fBins + fBins/2 + m];
			ijk.real[n*fBins + fBins/2 + m] = ijkShiftTemp.real[n*fBins + m];
			ijk.imag[n*fBins + fBins/2 + m] = ijkShiftTemp.imag[n*fBins + m];
		}
	}
	//Renormalize
	CVector normalization = CVector(s,1);
	normalization.real[0] = 1.0/nElements;
	ijk = ijk*normalization;
	return ijk;
}

CVector CVector::ijk2ijw(){ 
	CVector ijw = CVector(s);
	for (int n = 0; n < xBins; n++){
		for (int m = 0; m < yBins; m++){
			//Pull out the particular line of sight
			CVector f(s,fBins);
			for (int l = 0; l < fBins; l++){
				f.real[l] = real[l + fBins*m + fBins*yBins*n];
				f.imag[l] = imag[l + fBins*m + fBins*yBins*n];
			}
			//Shift
			CVector fShift(s,fBins);
			for (int l = 0; l < fBins; l++){
				fShift.real[l] = f.real[(fBins/2 + l)%fBins];
				fShift.imag[l] = f.imag[(fBins/2 + l)%fBins];
			}
			//Perform the FFT
			fftw_complex *in, *out;
			in = (fftw_complex*) fftw_malloc(fBins * sizeof(fftw_complex));
			out = (fftw_complex*) fftw_malloc(fBins * sizeof(fftw_complex));
			for (int l = 0; l < fBins; l++){
				in[l][0] = fShift.real[l];
				in[l][1] = fShift.imag[l];
			}
			fftw_plan FT = fftw_plan_dft_1d(fBins, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_execute(FT);
			CVector wShift(s,fBins);
			for (int l = 0; l < fBins; l++){
				wShift.real[l] = out[l][0];
				wShift.imag[l] = out[l][1];
			}
			fftw_destroy_plan(FT);
			fftw_free(in); fftw_free(out);
			//Shift
			CVector w(s,fBins);
			for (int l = 0; l < fBins; l++){
				w.real[l] = wShift.real[(fBins/2 + l)%fBins];
				w.imag[l] = wShift.imag[(fBins/2 + l)%fBins];
			}
			//Put back into original vector
			for (int l = 0; l < fBins; l++){
				ijw.real[l + fBins*m + fBins*yBins*n] = w.real[l];
				ijw.imag[l + fBins*m + fBins*yBins*n] = w.imag[l];
			}
		}
	}
	return ijw;
}

CVector CVector::ijw2ijk(){
		
	CVector ijk = CVector(s);
	for (int n = 0; n < xBins; n++){
		for (int m = 0; m < yBins; m++){
			//Pull out the particular line of sight
			CVector w(s,fBins);
			
			for (int l = 0; l < fBins; l++){
				w.real[l] = real[l + fBins*m + fBins*yBins*n];
				w.imag[l] = imag[l + fBins*m + fBins*yBins*n];
			}
			
			//Shift to start with the origin
			CVector wShift(s,fBins);
			for (int l = 0; l < fBins; l++){
				wShift.real[l] = w.real[(fBins/2 + l)%fBins];
				wShift.imag[l] = w.imag[(fBins/2 + l)%fBins];
			}
			
			//Perform the FFT
			fftw_complex *in, *out;
			in = (fftw_complex*) fftw_malloc(fBins * sizeof(fftw_complex));
			out = (fftw_complex*) fftw_malloc(fBins * sizeof(fftw_complex));
			for (int l = 0; l < fBins; l++){
				in[l][0] = wShift.real[l];
				in[l][1] = wShift.imag[l];
			}
			fftw_plan FT = fftw_plan_dft_1d(fBins, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
			fftw_execute(FT);
			CVector fShift(s,fBins);
			for (int l = 0; l < fBins; l++){
				fShift.real[l] = out[l][0];
				fShift.imag[l] = out[l][1];
			}
			fftw_destroy_plan(FT);
			fftw_free(in); fftw_free(out);
			//Shift
			CVector f(s,fBins);
			for (int l = 0; l < fBins; l++){
				f.real[l] = fShift.real[(fBins/2 + l)%fBins];
				f.imag[l] = fShift.imag[(fBins/2 + l)%fBins];
			}
			//Normalize and put back into original vector
			CVector normalization = CVector(s,1);
			normalization.real[0] = 1.0/fBins;
			f = f*normalization;	
			for (int l = 0; l < fBins; l++){
				ijk.real[l + fBins*m + fBins*yBins*n] = f.real[l];
				ijk.imag[l + fBins*m + fBins*yBins*n] = f.imag[l];
			}
		}
	}

	return ijk;

}

CVector CVector::ijk2uvk(){
	CVector uvk = CVector(s);
	for (int n = 0; n < fBins; n++){
		CVector ij(s,xBins*yBins);
		//Pull out a particular frequency slice
		for (int m = 0; m < xBins; m++){
			for (int l = 0; l < yBins; l++){
				ij.real[m*yBins + l] = real[n + l*fBins + m*fBins*yBins];
				ij.imag[m*yBins + l] = imag[n + l*fBins + m*fBins*yBins];				
			}
		}
		//Shift to start with origin
		CVector ijShift(s,xBins*yBins);
		for (int m = 0; m < xBins; m++){
			for (int l = 0; l < yBins; l++){
				ijShift.real[m*yBins + l] = ij.real[(yBins/2 + l)%yBins + m*yBins];
				ijShift.imag[m*yBins + l] = ij.imag[(yBins/2 + l)%yBins + m*yBins];
			}
		}
		ij = ijShift;
		for (int m = 0; m < xBins*yBins; m++){
			ijShift.real[m] = ij.real[(xBins*yBins/2 + m)%(xBins*yBins)];
			ijShift.imag[m] = ij.imag[(xBins*yBins/2 + m)%(xBins*yBins)];
		
		}
		//Fourier transform	
		fftw_complex *in, *out;
		in = (fftw_complex*) fftw_malloc(xBins*yBins * sizeof(fftw_complex));
		out = (fftw_complex*) fftw_malloc(xBins*yBins * sizeof(fftw_complex));
		for (int m = 0; m < xBins*yBins; m++){
			in[m][0] = ijShift.real[m];
			in[m][1] = ijShift.imag[m];
		}
		fftw_plan FT = fftw_plan_dft_2d(xBins, yBins, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(FT);
		CVector uvShift(s,xBins*yBins);
		for (int m = 0; m < xBins*yBins; m++){
			uvShift.real[m] = out[m][0];
			uvShift.imag[m] = out[m][1];
		}
		fftw_destroy_plan(FT);
		fftw_free(in); fftw_free(out);
		//Shift back
		CVector uv(s,xBins*yBins);
		for (int m = 0; m < xBins; m++){
			for (int l = 0; l < yBins; l++){
				uv.real[m*yBins + l] = uvShift.real[(yBins/2 + l)%yBins + m*yBins];
				uv.imag[m*yBins + l] = uvShift.imag[(yBins/2 + l)%yBins + m*yBins];
			}
		}
		uvShift = uv;
		for (int m = 0; m < xBins*yBins; m++){
			uv.real[m] = uvShift.real[(xBins*yBins/2 + m)%(xBins*yBins)];
			uv.imag[m] = uvShift.imag[(xBins*yBins/2 + m)%(xBins*yBins)];
		
		}
		//Put back into original vector
		for (int m = 0; m < xBins; m++){
			for (int l = 0; l < yBins; l++){
				uvk.real[n + l*fBins + m*fBins*yBins] = uv.real[m*yBins + l];
				uvk.imag[n + l*fBins + m*fBins*yBins] = uv.imag[m*yBins + l];
			}
		}
	}
	return uvk;
}

CVector CVector::uvk2ijk(){
	CVector ijk = CVector(s);
	for (int n = 0; n < fBins; n++){
		CVector uv(s,xBins*yBins);
		//Pull out a particular frequency slice
		for (int m = 0; m < xBins; m++){
			for (int l = 0; l < yBins; l++){
				uv.real[m*yBins + l] = real[n + l*fBins + m*fBins*yBins];
				uv.imag[m*yBins + l] = imag[n + l*fBins + m*fBins*yBins];				
			}
		}
		//Shift to start with origin
		CVector uvShift(s,xBins*yBins);
		for (int m = 0; m < xBins; m++){
			for (int l = 0; l < yBins; l++){
				uvShift.real[m*yBins + l] = uv.real[(yBins/2 + l)%yBins + m*yBins];
				uvShift.imag[m*yBins + l] = uv.imag[(yBins/2 + l)%yBins + m*yBins];
			}
		}
		uv = uvShift;
		for (int m = 0; m < xBins*yBins; m++){
			uvShift.real[m] = uv.real[(xBins*yBins/2 + m)%(xBins*yBins)];
			uvShift.imag[m] = uv.imag[(xBins*yBins/2 + m)%(xBins*yBins)];
		
		}
		//Fourier transform	
		fftw_complex *in, *out;
		in = (fftw_complex*) fftw_malloc(xBins*yBins * sizeof(fftw_complex));
		out = (fftw_complex*) fftw_malloc(xBins*yBins * sizeof(fftw_complex));
		for (int m = 0; m < xBins*yBins; m++){
			in[m][0] = uvShift.real[m];
			in[m][1] = uvShift.imag[m];
		}
		fftw_plan FT = fftw_plan_dft_2d(xBins, yBins, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(FT);
		CVector ijShift(s,xBins*yBins);
		for (int m = 0; m < xBins*yBins; m++){
			ijShift.real[m] = out[m][0];
			ijShift.imag[m] = out[m][1];
		}
		fftw_destroy_plan(FT);
		fftw_free(in); fftw_free(out);
		//Shift back
		CVector ij(s,xBins*yBins);
		for (int m = 0; m < xBins; m++){
			for (int l = 0; l < yBins; l++){
				ij.real[m*yBins + l] = ijShift.real[(yBins/2 + l)%yBins + m*yBins];
				ij.imag[m*yBins + l] = ijShift.imag[(yBins/2 + l)%yBins + m*yBins];
			}
		}
		ijShift = ij;
		for (int m = 0; m < xBins*yBins; m++){
			ij.real[m] = ijShift.real[(xBins*yBins/2 + m)%(xBins*yBins)];
			ij.imag[m] = ijShift.imag[(xBins*yBins/2 + m)%(xBins*yBins)];
		
		}
		//Normalize and put back into original vector
		CVector normalization = CVector(s,1);
		normalization.real[0] = 1.0/(xBins*yBins);
		ij = ij*normalization;	
		for (int m = 0; m < xBins; m++){
			for (int l = 0; l < yBins; l++){
				ijk.real[n + l*fBins + m*fBins*yBins] = ij.real[m*yBins + l];
				ijk.imag[n + l*fBins + m*fBins*yBins] = ij.imag[m*yBins + l];
			}
		}
	}
	return ijk;
}

CVector CVector::ij2uv(){
	//Shift to start with origin
	CVector ijShift(s,xBins*yBins);
	for (int m = 0; m < xBins; m++){
		for (int l = 0; l < yBins; l++){
			ijShift.real[m*yBins + l] = real[(yBins/2 + l)%yBins + m*yBins];
			ijShift.imag[m*yBins + l] = imag[(yBins/2 + l)%yBins + m*yBins];
		}
	}
	CVector ij = ijShift;
	for (int m = 0; m < xBins*yBins; m++){
		ijShift.real[m] = ij.real[(xBins*yBins/2 + m)%(xBins*yBins)];
		ijShift.imag[m] = ij.imag[(xBins*yBins/2 + m)%(xBins*yBins)];
	
	}
	//Fourier transform	
	fftw_complex *in, *out;
	in = (fftw_complex*) fftw_malloc(xBins*yBins * sizeof(fftw_complex));
	out = (fftw_complex*) fftw_malloc(xBins*yBins * sizeof(fftw_complex));
	for (int m = 0; m < xBins*yBins; m++){
		in[m][0] = ijShift.real[m];
		in[m][1] = ijShift.imag[m];
	}
	fftw_plan FT = fftw_plan_dft_2d(xBins, yBins, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(FT);
	CVector uvShift(s,xBins*yBins);
	for (int m = 0; m < xBins*yBins; m++){
		uvShift.real[m] = out[m][0];
		uvShift.imag[m] = out[m][1];
	}
	fftw_destroy_plan(FT);
	fftw_free(in); fftw_free(out);
	//Shift back
	CVector uv(s,xBins*yBins);
	for (int m = 0; m < xBins; m++){
		for (int l = 0; l < yBins; l++){
			uv.real[m*yBins + l] = uvShift.real[(yBins/2 + l)%yBins + m*yBins];
			uv.imag[m*yBins + l] = uvShift.imag[(yBins/2 + l)%yBins + m*yBins];
		}
	}
	uvShift = uv;
	for (int m = 0; m < xBins*yBins; m++){
		uv.real[m] = uvShift.real[(xBins*yBins/2 + m)%(xBins*yBins)];
		uv.imag[m] = uvShift.imag[(xBins*yBins/2 + m)%(xBins*yBins)];
	
	}
	return uv;
}

CVector CVector::uv2ij(){
	//Shift to start with origin
	CVector uvShift(s,xBins*yBins);
	for (int m = 0; m < xBins; m++){
		for (int l = 0; l < yBins; l++){
			uvShift.real[m*yBins + l] = real[(yBins/2 + l)%yBins + m*yBins];
			uvShift.imag[m*yBins + l] = imag[(yBins/2 + l)%yBins + m*yBins];
		}
	}
	CVector uv = uvShift;
	for (int m = 0; m < xBins*yBins; m++){
		uvShift.real[m] = uv.real[(xBins*yBins/2 + m)%(xBins*yBins)];
		uvShift.imag[m] = uv.imag[(xBins*yBins/2 + m)%(xBins*yBins)];
	
	}
	//Fourier transform	
	fftw_complex *in, *out;
	in = (fftw_complex*) fftw_malloc(xBins*yBins * sizeof(fftw_complex));
	out = (fftw_complex*) fftw_malloc(xBins*yBins * sizeof(fftw_complex));
	for (int m = 0; m < xBins*yBins; m++){
		in[m][0] = uvShift.real[m];
		in[m][1] = uvShift.imag[m];
	}
	fftw_plan FT = fftw_plan_dft_2d(xBins, yBins, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(FT);
	CVector ijShift(s,xBins*yBins);
	for (int m = 0; m < xBins*yBins; m++){
		ijShift.real[m] = out[m][0];
		ijShift.imag[m] = out[m][1];
	}
	fftw_destroy_plan(FT);
	fftw_free(in); fftw_free(out);
	//Shift back
	CVector ij(s,xBins*yBins);
	for (int m = 0; m < xBins; m++){
		for (int l = 0; l < yBins; l++){
			ij.real[m*yBins + l] = ijShift.real[(yBins/2 + l)%yBins + m*yBins];
			ij.imag[m*yBins + l] = ijShift.imag[(yBins/2 + l)%yBins + m*yBins];
		}
	}
	ijShift = ij;
	for (int m = 0; m < xBins*yBins; m++){
		ij.real[m] = ijShift.real[(xBins*yBins/2 + m)%(xBins*yBins)]/(xBins*yBins); //Normalize
		ij.imag[m] = ijShift.imag[(xBins*yBins/2 + m)%(xBins*yBins)]/(xBins*yBins);
	
	}
	
	return ij;
}


CVector CVector::zeroPadThis(){
	CVector result(s,xBins*yBins*fBins*zeroPad*zeroPad*zeroPad);
	for (int i = 0; i < xBins; i++){
		int xThere = i + (zeroPad-1)*xBins/2;
		for (int j = 0; j < yBins; j++){
			int yThere = j + (zeroPad-1)*yBins/2;
			for (int k = 0; k < fBins; k++){
				int fThere = k + (zeroPad-1)*fBins/2;
				int there = xThere*yBins*zeroPad*fBins*zeroPad + yThere*fBins*zeroPad + fThere;
				int here = i*yBins*fBins + j*fBins + k;
				result.real[there] = real[here];
				result.imag[there] = imag[here];
			}
		}
	}
	
	return result;
}



CVector CVector::FTforPowerSpectrum(bool forward){
	xBins *= zeroPad;
	yBins *= zeroPad;
	fBins *= zeroPad;

	//Shift to start with origin
	CVector shiftedVector(s,nElements);	
	if (nElements != xBins*yBins*fBins){
		cout << endl << "WARNING: Vector to be FTTed for Power Spectrum is the wrong size!" << endl << endl;
		return shiftedVector;
	}	
	for (int i = 0; i < xBins; i++){
		int xNew = (i + xBins/2)%xBins; //"New" means the final location for a given initial location (shifted over)
		for (int j = 0; j < yBins; j++){
			int yNew = (j + yBins/2)%yBins;
			for (int k = 0; k < fBins; k++){
				int fNew = (k + fBins/2)%fBins;
				shiftedVector.real[xNew*yBins*fBins + yNew*fBins + fNew] = real[i*yBins*fBins + j*fBins + k];
				shiftedVector.imag[xNew*yBins*fBins + yNew*fBins + fNew] = imag[i*yBins*fBins + j*fBins + k];
			}
		}
	}
	
	//Perform the FFT (this is all FFTW boilerplate)
	fftw_complex *in, *out;
	in = (fftw_complex*) fftw_malloc(xBins*yBins*fBins * sizeof(fftw_complex));
	out = (fftw_complex*) fftw_malloc(xBins*yBins*fBins * sizeof(fftw_complex));
	for (int n = 0; n < shiftedVector.nElements; n++){
		in[n][0] = shiftedVector.real[n];
		in[n][1] = shiftedVector.imag[n];
	}
	fftw_plan fourierTransform;
	if (forward) fourierTransform = fftw_plan_dft_3d(xBins, yBins, fBins, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	else fourierTransform = fftw_plan_dft_3d(xBins, yBins, fBins, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(fourierTransform);
	CVector FTShift(s,nElements);
	for (int n = 0; n < FTShift.nElements; n++){
		FTShift.real[n] = out[n][0];
		FTShift.imag[n] = out[n][1];
	}
	fftw_destroy_plan(fourierTransform);
	fftw_free(in); fftw_free(out);
	
	//Shift back (exact same procedure as before)
	CVector FT(s,nElements);
	for (int i = 0; i < xBins; i++){
		int xNew = (i + xBins/2)%xBins;
		for (int j = 0; j < yBins; j++){
			int yNew = (j + yBins/2)%yBins;
			for (int k = 0; k < fBins; k++){
				int fNew = (k + fBins/2)%fBins;
				FT.real[xNew*yBins*fBins + yNew*fBins + fNew] = FTShift.real[i*yBins*fBins + j*fBins + k];
				FT.imag[xNew*yBins*fBins + yNew*fBins + fNew] = FTShift.imag[i*yBins*fBins + j*fBins + k];
			}
		}
	}
	xBins /= zeroPad;
	yBins /= zeroPad;
	fBins /= zeroPad;
	
	return FT;
}

CVector CVector::zeroPadAndFTforPowerSpectrum(bool forward){
	CVector vector1(s,xBins*yBins*fBins*zeroPad*zeroPad*zeroPad);
	CVector vector2(s,xBins*yBins*fBins*zeroPad*zeroPad*zeroPad);
	
	for (int i = 0; i < xBins; i++){
		int xThere = i + (zeroPad-1)*xBins/2;
		for (int j = 0; j < yBins; j++){
			int yThere = j + (zeroPad-1)*yBins/2;
			for (int k = 0; k < fBins; k++){
				int fThere = k + (zeroPad-1)*fBins/2;
				int there = xThere*yBins*zeroPad*fBins*zeroPad + yThere*fBins*zeroPad + fThere;
				int here = i*yBins*fBins + j*fBins + k;
				vector1.real[there] = real[here]; //zeroPadded vector
				vector1.imag[there] = imag[here];
			}
		}
	}
	
	xBins *= zeroPad;
	yBins *= zeroPad;
	fBins *= zeroPad;

	//Shift to start with origin
	for (int i = 0; i < xBins; i++){
		int xNew = (i + xBins/2)%xBins; //"New" means the final location for a given initial location (shifted over)
		for (int j = 0; j < yBins; j++){
			int yNew = (j + yBins/2)%yBins;
			for (int k = 0; k < fBins; k++){
				int fNew = (k + fBins/2)%fBins;
				vector2.real[xNew*yBins*fBins + yNew*fBins + fNew] = vector1.real[i*yBins*fBins + j*fBins + k]; //shifted and zeroPadded vector
				vector2.imag[xNew*yBins*fBins + yNew*fBins + fNew] = vector1.imag[i*yBins*fBins + j*fBins + k];
			}
		}
	}
	
	
	//Perform the FFT (this is all FFTW boilerplate)
	fftw_complex *FTarray;
	FTarray = (fftw_complex*) fftw_malloc(xBins*yBins*fBins * sizeof(fftw_complex));
	for (int n = 0; n < xBins*yBins*fBins; n++){
		FTarray[n][0] = vector2.real[n];
		FTarray[n][1] = vector2.imag[n];
	}
	
	fftw_plan fourierTransform;
	if (forward) fourierTransform = fftw_plan_dft_3d(xBins, yBins, fBins, FTarray, FTarray, FFTW_FORWARD, FFTW_ESTIMATE);
	else fourierTransform = fftw_plan_dft_3d(xBins, yBins, fBins, FTarray, FTarray, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(fourierTransform);
	for (int n = 0; n < xBins*yBins*fBins; n++){
		vector1.real[n] = FTarray[n][0]; //shifted, FT'ed, and ZP'ed array
		vector1.imag[n] = FTarray[n][1];
	}
	fftw_destroy_plan(fourierTransform);
	fftw_free(FTarray);
	
	//Shift back (exact same procedure as before)
	for (int i = 0; i < xBins; i++){
		int xNew = (i + xBins/2)%xBins;
		for (int j = 0; j < yBins; j++){
			int yNew = (j + yBins/2)%yBins;
			for (int k = 0; k < fBins; k++){
				int fNew = (k + fBins/2)%fBins;
				vector2.real[xNew*yBins*fBins + yNew*fBins + fNew] = vector1.real[i*yBins*fBins + j*fBins + k]; //FT'ed, and ZP'ed array
				vector2.imag[xNew*yBins*fBins + yNew*fBins + fNew] = vector1.imag[i*yBins*fBins + j*fBins + k];
			}
		}
	}
	xBins /= zeroPad;
	yBins /= zeroPad;
	fBins /= zeroPad;
	
	return vector2;
}

//Reorders from rising with frequency to falling with frequency or vice versa
CVector CVector::reverseFrequencyComponents(){
	CVector reversed = CVector(s);
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			for (int k = 0; k < fBins; k++){
				reversed.real[i*yBins*fBins + j*fBins + k] = real[i*yBins*fBins + j*fBins + (fBins - k - 1)];
				reversed.imag[i*yBins*fBins + j*fBins + k] = imag[i*yBins*fBins + j*fBins + (fBins - k - 1)];
			}
		}
	}	
	return reversed;
}

