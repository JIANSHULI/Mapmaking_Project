#include "Specs.hpp"
#ifndef COMPLEX_VECTOR
#define COMPLEX_VECTOR

using namespace std;

class CVector {

    public:
		CVector();
		CVector(Specs *s);
		CVector(Specs *s, int N);
		CVector(Specs *s, string filename);
		CVector(const CVector& copy);
        ~CVector();
		void loadSpecs(Specs *s);
		CVector dot(const CVector& x);
		CVector& operator=(const CVector& x);
		CVector operator*(const CVector& c);
		CVector operator+(const CVector& x);
		CVector operator-(const CVector& x);
		CVector operator/(const CVector& d);
		bool allMagnitudesLessThan(double c);
		double magnitude();
        void printAll(string name);
		void printRealToFile(string filename);
		CVector ijk2uvw();
		CVector uvw2ijk();
		CVector ijk2ijw();
		CVector ijw2ijk();
		CVector ijk2uvk();
		CVector uvk2ijk();
		CVector ij2uv();
		CVector uv2ij();
		CVector zeroPadThis();
		CVector FTforPowerSpectrum(bool forward);
		CVector zeroPadAndFTforPowerSpectrum(bool forward);
		CVector reverseFrequencyComponents();
        int nElements;
		int xBins;
		int yBins;
		int fBins;
		int zeroPad;
        double *real;
        double *imag;
		Specs *s;
};

#endif

