#include "CVector.hpp"
#include "Specs.hpp"
#include <vector>
using namespace std;

class Toeplitz {

    public:
		Toeplitz();
		Toeplitz(int n, string filename, string dim);
		Toeplitz(int n, string filename, int x, int y);
		Toeplitz(const Toeplitz& copy);
		Toeplitz& operator=(const Toeplitz& x);
        ~Toeplitz();
		Toeplitz operator+(const Toeplitz& T);
        CVector operator*(const CVector& x);
		vector< vector<double> > gaussRandField2D(Toeplitz Y);
        void printAll();
        int nElements;
		bool xDim;
		bool yDim;
		bool fDim;
		int xPos;
		int yPos;
        double *entry;
		Specs *s;
		
};
