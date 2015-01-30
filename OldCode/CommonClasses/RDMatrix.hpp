#include "CVector.hpp"
#include "Specs.hpp"
using namespace std;

class RDMatrix {

    public:
        RDMatrix();
		RDMatrix(Specs *s);
		RDMatrix(Specs *s, int N);
		RDMatrix(Specs *s, string filename);
		RDMatrix(const RDMatrix& copy);
		RDMatrix& operator=(const RDMatrix& x);
        ~RDMatrix();
		void loadSpecs(Specs *s);
        CVector operator*(const CVector& x);
        void printAll();
        int nElements;
		int xBins;
		int yBins;
		int fBins;
		double rms;
        double *entry;
		Specs *s;

};
