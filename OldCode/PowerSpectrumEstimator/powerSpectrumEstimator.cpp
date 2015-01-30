#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "../CommonClasses/RDMatrix.hpp"
#include "../CommonClasses/CVector.hpp"
#include "../CommonClasses/Specs.hpp"
#include <stdio.h>
#include <time.h>
#include <vector>
#include "../CommonClasses/Toeplitz.hpp"
#include <string>
#include <cstring>
#include <sstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
#include "fftw3.h"
#include <stdlib.h>
#include <ctype.h>

using namespace std;

/*******************************************************************
CONSTANTS, GLOBAL VARIABLES, AND DATA STRUCTURES
********************************************************************/

const double pi = 3.1415926535897932384626433832795;
const double c = 299792000; //m/s
const double H0 = 67770; //m/s/Mpc
const double OmegaM = .3086;
const double OmegaL = .6914;
const double f21cm = 1420.41; //MHz
Specs *s;
int xBins, yBins, fBins, nElements;
double xyLength, fLength, fStart;
int kParaBins, kPerpBins, kSphereBins, zeroPad;
double CGBound, PreconEVThreshold, horizonWedgeBuffer;
int FisherMCNum;
bool verbose, turnOffPreconditioner, fullGammaEigenDecomp, removeMean, JEquals1, includeCroppedPerpBins, useCrossPower, projectOutHorizonWedge;
bool CequalsI, CequalsN, CequalsNU, removeR, removeG;
bool printDataCubesToFile, printAllPSE, saveKBinEdges;
bool MatrixRecon, QMatrixRecon, CovTesting;
int covTestingTrials, binLogarithmically;
int fields, centerField;
string qHatsDirectory = "../QuadraticEstimators/";
string sphericalQHatsDirectory = "../QuadraticEstimators/Spherical/";
string kParaBinCentersFilename = "../QuadraticEstimators/kParaBinCenters.dat";
string kPerpBinCentersFilename = "../QuadraticEstimators/kPerpBinCenters.dat";
string kSphericalBinCentersFilename = "../QuadraticEstimators/Spherical/kSphericalBinCenters.dat";
string kParaBinEdgesFilename = "../QuadraticEstimators/kParaBinEdges.dat";
string kPerpBinEdgesFilename = "../QuadraticEstimators/kPerpBinEdges.dat";
string kSphericalBinEdgesFilename = "../QuadraticEstimators/Spherical/kSphericalBinEdges.dat";
string cubeDirectory = "../Cubes/";
string dataCubeMaskFilename = "../Cubes/dataCubeMask.dat";
string dataCubeString1 = "_field_";
string datacubeString2 = "_slice_";
string datacubeString3 = ".dat";


struct ResolvedCovariance {
	int nSources;
	double spectralIndexPercentError, fluxPercentError, omegaPix, referenceFrequency;
	vector< vector<int> > positions;
	vector<double> fluxes, spectralIndices, freqs;
	vector<Toeplitz> Cov;
	vector< vector< vector<double> > > Eigen;
	vector< vector< vector<double> > > EigenMasked;
	RDMatrix Dcross, Dmean;	
	vector< vector<int> > eigenToUseInGamma;
};

struct UnresolvedCovariance {
	double spectralIndexMean, spectralIndexStd, referenceFrequency, omegaPix, I1, I2;
	vector<double> freqs;
	Toeplitz xCov, yCov, fCov;
	RDMatrix Dcross, Dmean;	
	vector< vector<double> > fEigen;
	vector< vector<double> > xEigen; 
	vector< vector<double> > yEigen;
	vector< vector<double> > fdEigen;
	vector< vector<double> > fEigenMasked;
	vector< vector<double> > fdEigenMasked;
};

struct GalacticCovariance {
	double spectralIndexMean, spectralIndexStd, referenceFrequency, omegaPix, I1, I2;
	vector<double> freqs;
	Toeplitz xCov, yCov, fCov;
	RDMatrix Dcross, Dmean;	
	vector< vector<double> > fEigen;
	vector< vector<double> > xEigen; 
	vector< vector<double> > yEigen; 
	vector< vector<double> > fdEigen;
	vector< vector<double> > fEigenMasked;
	vector< vector<double> > fdEigenMasked;
	vector< vector< vector<int> > > eigenToUseInGamma;
};

struct NoiseCovariance {
	RDMatrix Cov, effCov, avgCov;
};

struct Preconditioner {
	double noiseMin;
	vector< vector< CVector > > GammaEigenvectorsFourier; //storage: [k (may not go up to fBins)][eigenvalue/vector number (at least 1)].real/imag[eigenvector entry];
	vector< vector< CVector > > GammaEigenvectorsFourierTrans;
	vector< vector<double> > GammaEigenvalues; //storage: [k (may not go up to fBins)][eigenvalue]
	vector< vector<double> > GammaBarFactors;
	CVector mask;
};

/*******************************************************************
LOADING DATA
********************************************************************/

void loadSpecs(string dataSpecsFilename, string cubeParametersFilename, string powerSpectrumSpecsFilename, string powerSpectrumOptionsFilename){
	//Load in locations of real data cube and its parameters
	fstream infile(dataSpecsFilename.c_str(),fstream::in);
	string dummy;
	for (int n = 0; n < 4; n++) infile >> dummy >> dummy;
	infile >> dummy >> fields;
	infile >> dummy >> centerField;
	infile.close();
	
	//Load in data cube parameters
	infile.open(cubeParametersFilename.c_str(),fstream::in);
	infile >> dummy >> xBins;
	infile >> dummy >> yBins;
	infile >> dummy >> fBins;
	infile >> dummy >> xyLength;
	infile >> dummy >> fLength;
	infile >> dummy >> fStart; 
	infile.close();
	
	//Load in power spectrum estimation specs
	infile.open(powerSpectrumSpecsFilename.c_str(),fstream::in);
	infile >> dummy >> zeroPad;
	infile >> dummy >> PreconEVThreshold;
	infile >> dummy >> CGBound;
	infile >> dummy >> FisherMCNum;
	infile.close();
	
	//Load in advanced power spectrum estimation options
	infile.open(powerSpectrumOptionsFilename.c_str(),fstream::in);
	infile >> dummy >> dummy >> verbose;
	infile >> dummy >> dummy >> turnOffPreconditioner;
	infile >> dummy >> dummy >> fullGammaEigenDecomp;
	infile >> dummy >> dummy >> CequalsI;
	infile >> dummy >> dummy >> CequalsN;
	infile >> dummy >> dummy >> CequalsNU;
	infile >> dummy >> dummy >> removeR;
	infile >> dummy >> dummy >> removeG;
	infile >> dummy >> dummy >> removeMean;
	infile >> dummy >> dummy >> printDataCubesToFile;
	infile >> dummy >> dummy >> CovTesting;
	infile >> dummy >> dummy >> covTestingTrials;
	infile >> dummy >> dummy >> MatrixRecon;
	infile >> dummy >> dummy >> QMatrixRecon;
	infile >> dummy >> dummy >> printAllPSE;
	infile >> dummy >> dummy >> saveKBinEdges;
	infile >> dummy >> dummy >> JEquals1;
	infile >> dummy >> dummy >> includeCroppedPerpBins;
	infile >> dummy >> dummy >> useCrossPower;
	infile >> dummy >> dummy >> binLogarithmically;
	infile >> dummy >> dummy >> projectOutHorizonWedge;
	infile >> dummy >> dummy >> horizonWedgeBuffer;
	infile.close();
	
	//Load relevant specifications into Specs object
	s->xBins = xBins;
	s->yBins = yBins;
	s->fBins = fBins;
	s->zeroPad = zeroPad;
	s->CGBound = CGBound;
	s->xyLength = xyLength;
	s->fLength = fLength;
	s->fStart = fStart;
	s->PreconEVThreshold = PreconEVThreshold;
	s->FisherMCNum = FisherMCNum;
	nElements = xBins*yBins*fBins;
	kParaBins = fBins/2 + 1;
	kSphereBins = kParaBins + 1;
	if (includeCroppedPerpBins){
		kPerpBins = int(round(sqrt(pow(xBins/2,2)+pow(yBins/2,2)))) + 1;
	} else {
		kPerpBins = xBins/2 + 1;
	}

}

CVector loadRealData(string filename){
	CVector realData(s);
	fstream infile(filename.c_str(),fstream::in);
	double value = 0.0;
	for (int n = 0; n < nElements; n++){
		infile >> value;
		realData.real[n] = value;
	}
	return realData;
}

//Loads a file in my format into a CVector
CVector loadFileIntoCVector(int field, int slice, string type){
	stringstream ss;
	ss << cubeDirectory << type << dataCubeString1 << field << datacubeString2 << slice << datacubeString3;
	fstream infile((ss.str()).c_str(),fstream::in);
	CVector loadedVector(s);
	double value;
	for (int n = 0; n < nElements; n++){
		infile >> value;
		loadedVector.real[n] = value;
	}
	return loadedVector;
}

//Loads in the eigenvectors and eigenvalues of the Toeplitz matrices for gaussian random data
vector< vector<double> > loadToeplitzEigensystem(int size, string filename){
	fstream infile(filename.c_str(),fstream::in);
	int numEVs;
	infile >> numEVs;
	vector< vector<double> > eigensystem(numEVs, vector<double> (size + 1,0));
	for (int n = 0; n < numEVs; n++){
		for (int m = 0; m < size + 1; m++){
			infile >> eigensystem[n][m];
		}
	}
	return eigensystem;
}

//Loads all the small Toeplitz matrices for each resolved point source into a vector of such matrices
void loadRToeplitz(ResolvedCovariance& R, string filePrefix, string fileSuffix){
	for (int n = 0; n < R.nSources; n++){
		stringstream ss;
		ss << filePrefix << n << fileSuffix;	
		Toeplitz RT(fBins,ss.str(),R.positions[n][0],R.positions[n][1]);
		R.Cov.push_back(RT);
	}
}

//Loads in all the toeplitz eigensystems for each resolved point source
vector< vector< vector<double> > > loadRToeplitzEigensystems(ResolvedCovariance& R, int size, int nSources, string filePrefix, string fileSuffix, bool masked){
	vector< vector< vector<double> > > eigensystems(nSources, vector< vector<double> >(size, vector<double> (size + 1,0))); //TODO: This vector could be smaller, because it has room for as many EVs as possible
	for (int n = 0; n < nSources; n++){
		stringstream ss;
		ss << filePrefix << n << fileSuffix;
		string filename = ss.str();
		if (masked) {
			R.EigenMasked.push_back(loadToeplitzEigensystem(size,filename));
		} else {
			R.Eigen.push_back(loadToeplitzEigensystem(size,filename));
		}
	}
	return eigensystems;
}



void loadR(ResolvedCovariance& R, string RSourcesFile, string RToeplitzPrefix, string RToeplitzSuffix, string REigenPrefix , string REigenSuffix, string REigenMaskedPrefix, string REigenMaskedSuffix, string RDcrossFile, string RDmeanFile){	
	//Load in Source Info
	fstream infile(RSourcesFile.c_str(),fstream::in);
	infile >> R.nSources >> R.spectralIndexPercentError >> R.fluxPercentError >> R.referenceFrequency >> R.omegaPix;
	vector < vector<int> > positions(R.nSources, vector<int>(2,0));
	vector<double> fluxes(R.nSources,0);
	vector<double> spectralIndices(R.nSources,0);
	vector<double> freqs(fBins,0);
	for (int n = 0; n < R.nSources; n++){
		infile >> positions[n][0] >> positions[n][1] >> fluxes[n] >> spectralIndices[n];
	}
	for (int k = 0; k < fBins; k++) infile >> freqs[k];
	R.positions = positions;
	R.fluxes = fluxes;
	R.spectralIndices = spectralIndices;
	R.freqs = freqs;
	infile.close();

	//Load in diagonal Matrices
	R.Dcross = RDMatrix(s, RDcrossFile);
	R.Dmean = RDMatrix(s, RDmeanFile);

	//Load in Toeplitz Matrices and Eigensystems
	loadRToeplitz(R, RToeplitzPrefix, RToeplitzSuffix);
	loadRToeplitzEigensystems(R, fBins, R.nSources, REigenPrefix, REigenSuffix, false);
	loadRToeplitzEigensystems(R, fBins, R.nSources, REigenMaskedPrefix, REigenMaskedSuffix, true);
}


void loadU(UnresolvedCovariance& U, string UInfoFile, string UXToeplitzFile, string UYToeplitzFile, string UFToeplitzFile, string UDcrossFile, string UDmeanFile, string UXEigenFile, string UYEigenFile, string UFEigenFile, string UFDEigenFile, string UFEigenMaskedFile, string UFDEigenMaskedFile){
	//Load Unresolved Info
	fstream infile(UInfoFile.c_str(),fstream::in);
	infile >> U.spectralIndexMean >> U.spectralIndexStd >> U.referenceFrequency >> U.omegaPix >> U.I1 >> U.I2;
	vector<double> freqs(fBins,0);
	for (int k = 0; k < fBins; k++) infile >> freqs[k];
	U.freqs = freqs;
	infile.close();
	
	//Load in Toeplitz Matrices
	U.xCov = Toeplitz(xBins, UXToeplitzFile, "x");
	U.yCov = Toeplitz(yBins, UYToeplitzFile, "y");
	U.fCov = Toeplitz(fBins, UFToeplitzFile, "f");
	
	//Load in Eigensystems
	U.xEigen = loadToeplitzEigensystem(xBins, UXEigenFile);
	U.yEigen = loadToeplitzEigensystem(yBins, UYEigenFile);
	U.fEigen = loadToeplitzEigensystem(fBins, UFEigenFile);
	U.fdEigen = loadToeplitzEigensystem(fBins, UFDEigenFile);
	U.fEigenMasked = loadToeplitzEigensystem(fBins, UFEigenMaskedFile);
	U.fdEigenMasked = loadToeplitzEigensystem(fBins, UFDEigenMaskedFile);
	
	//Load in Diagonal Matrices
	U.Dcross = RDMatrix(s,UDcrossFile);
	U.Dmean = RDMatrix(s, UDmeanFile);
}

void loadG(GalacticCovariance& G, string GInfoFile, string GXToeplitzFile, string GYToeplitzFile, string GFToeplitzFile, string GDcrossFile, string GDmeanFile, string GXEigenFile, string GYEigenFile, string GFEigenFile, string GFDEigenFile, string GFEigenMaskedFile, string GFDEigenMaskedFile){
	//Load Unresolved Info
	fstream infile(GInfoFile.c_str(),fstream::in);
	infile >> G.spectralIndexMean >> G.spectralIndexStd >> G.referenceFrequency >> G.omegaPix >> G.I1 >> G.I2;
	vector<double> freqs(fBins,0);
	for (int k = 0; k < fBins; k++) infile >> freqs[k];
	G.freqs = freqs;
	infile.close();
	
	//Load in Toeplitz Matrices
	G.xCov = Toeplitz(xBins, GXToeplitzFile, "x");
	G.yCov = Toeplitz(yBins, GYToeplitzFile, "y");
	G.fCov = Toeplitz(fBins, GFToeplitzFile, "f");
	
	//Load in Eigensystems
	G.xEigen = loadToeplitzEigensystem(xBins, GXEigenFile);
	G.yEigen = loadToeplitzEigensystem(yBins, GYEigenFile);
	G.fEigen = loadToeplitzEigensystem(fBins, GFEigenFile);
	G.fdEigen = loadToeplitzEigensystem(fBins, GFDEigenFile);
	G.fEigenMasked = loadToeplitzEigensystem(fBins, GFEigenMaskedFile);
	G.fdEigenMasked = loadToeplitzEigensystem(fBins, GFDEigenMaskedFile);
	
	//Load in Diagonal Matrices
	G.Dcross = RDMatrix(s, GDcrossFile);
	G.Dmean = RDMatrix(s, GDmeanFile);
}

void loadN(NoiseCovariance& N, string NFile, CVector mask){
	//Load in N 
	N.Cov = RDMatrix(s,NFile);
	//Calculate Neff
 	N.effCov = N.Cov;
	for (int n = 0; n < N.Cov.nElements; n++){
		if (N.Cov.entry[n] == -1 || mask.real[n] == 1) N.effCov.entry[n] = 1;  //These are the modes that will be projected out by the psuedoinverse
	}
	N.avgCov = N.effCov;
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			int count = 0;
			double sum = 0;
			for (int k = 0; k < fBins; k++){
				if (N.Cov.entry[i*fBins*yBins + j*fBins + k] > 0){
					count++;
					sum += N.effCov.entry[i*fBins*yBins + j*fBins + k];
				}
			}
			for (int k = 0; k < fBins; k++){
				if (count == 0) N.avgCov.entry[i*fBins*yBins + j*fBins + k] = 1;
				else if (mask.real[i*fBins*yBins + j*fBins + k] == 0) N.avgCov.entry[i*fBins*yBins + j*fBins + k] = sum/count;
			}
		}
	}

	//N.effCov = N.avgCov; //This should be commented out...it's only for testing to see if it makes preconditioning faster

/*	string outFilename1 = "NeffCov.dat";
	string outFilename2 = "NavgCov.dat";
	ofstream outfile1, outfile2;
	outfile1.precision(30);
	outfile2.precision(30);
	outfile1.open(outFilename1.c_str(), ios::trunc);
	outfile2.open(outFilename2.c_str(), ios::trunc);
	for (int n = 0; n < nElements; n++){
		outfile1 << N.effCov.entry[n] << endl;
		outfile2 << N.avgCov.entry[n] << endl;
	}
	outfile1.close();
	outfile2.close();*/

	return;
}

/*******************************************************************
CONSTRUCTING AND APPLYING THE PRECONDITIONER
********************************************************************/


void applyMask(CVector& dataCube, CVector& mask){
	for (int n = 0; n < nElements; n++){
		dataCube.real[n] *= (1.0 - mask.real[n]);
		dataCube.imag[n] *= (1.0 - mask.real[n]);
	}
}

//Does the full eigendecomposition of Gamma.  This might end up being faster and more accurate than the Lanczos algorithm.
int GammaEigenDecomp(int k, Preconditioner& precon, ResolvedCovariance& R, GalacticCovariance& G, NoiseCovariance& N){
	
	double* Gamma;
	Gamma= new double[xBins*yBins*xBins*yBins];
	for (int n = 0; n < xBins*yBins*xBins*yBins; n++){
		Gamma[n] = 0;
	}
	
	//Add Rxy into Gamma;
	if (!removeR){
		for (int m = 0; m < R.eigenToUseInGamma[k].size(); m++){
			int n = R.eigenToUseInGamma[k][m];
			if (R.eigenToUseInGamma[k][0] != -1){
				int here = R.positions[n][1] + R.positions[n][0]*yBins;
				double projection = 0;
				for (int k2 = 0; k2 < fBins; k2++) projection += R.EigenMasked[n][k][k2+1] * G.fdEigen[k][k2+1];
				if (projection < 0) projection = -projection;
				double lambdaR = projection*R.EigenMasked[n][k][0];

				CVector EigenVectorWithUnobservedBaselines(s,xBins*yBins);
				EigenVectorWithUnobservedBaselines.fBins = 1;
				EigenVectorWithUnobservedBaselines.real[here] = 1;

				CVector EigenVectorWithUnobservedBaselinesFT = EigenVectorWithUnobservedBaselines.ij2uv();
				int NCovChannelToUse = 0;
				while (precon.mask.real[NCovChannelToUse] == 1) NCovChannelToUse++;
				for (int i = 0; i < xBins; i++){
					for (int j = 0; j < yBins; j++){
						if (N.Cov.entry[i*yBins*fBins + j*fBins + NCovChannelToUse] == -1){
							EigenVectorWithUnobservedBaselinesFT.real[i*yBins + j] = 0;
							EigenVectorWithUnobservedBaselinesFT.imag[i*yBins + j] = 0;
						}
					}
				}
				CVector EigenVectorWithoutUnobservedBaselines = EigenVectorWithUnobservedBaselinesFT.uv2ij();

				for (int i = 0; i < xBins; i++){
					for (int j = 0; j < yBins; j++){
						for (int i2 = 0; i2 < xBins; i2++){
							for (int j2 = 0; j2 < yBins; j2++){
								int here1 = i*yBins + j;
								int here2 = i2*yBins + j2;
								int here = here1*xBins*yBins + here2;
								Gamma[here] += lambdaR * EigenVectorWithoutUnobservedBaselines.real[i*yBins + j] * EigenVectorWithoutUnobservedBaselines.real[i2*yBins + j2];
							}
						}
					}
				}
			}
		}
	} 


	//Add Gxy into Gamma
	if (!removeG){
		for (int m = 0; m < G.eigenToUseInGamma[k].size(); m++){
			if (G.eigenToUseInGamma[k][0][0] != -1){
				int xEig = G.eigenToUseInGamma[k][m][0];
				int yEig = G.eigenToUseInGamma[k][m][1];
				double lambdaG = G.xEigen[xEig][0]*G.yEigen[yEig][0]*G.fdEigenMasked[k][0];
				CVector EigenVectorWithUnobservedBaslines(s,xBins*yBins);
				EigenVectorWithUnobservedBaslines.fBins = 1;
				
				for (int i = 0; i < xBins; i++){
					for (int j = 0; j < yBins; j++){
						EigenVectorWithUnobservedBaslines.real[i*yBins + j] = G.xEigen[xEig][i+1] * G.yEigen[yEig][j+1];
					}
				}

				
				CVector EigenVectorWithUnobservedBaslinesFT = EigenVectorWithUnobservedBaslines.ij2uv();
				
				int NCovChannelToUse = 0;
				while (precon.mask.real[NCovChannelToUse] == 1) NCovChannelToUse++;
				for (int i = 0; i < xBins; i++){
					for (int j = 0; j < yBins; j++){
						if (N.Cov.entry[i*yBins*fBins + j*fBins + NCovChannelToUse] == -1){
							EigenVectorWithUnobservedBaslinesFT.real[i*yBins + j] = 0;
							EigenVectorWithUnobservedBaslinesFT.imag[i*yBins + j] = 0;
						}
					}
				}
				CVector EigenVectorWithoutUnobservedBaslines = 	EigenVectorWithUnobservedBaslinesFT.uv2ij();
				
				//(EigenVectorWithoutUnobservedBaslines-EigenVectorWithUnobservedBaslines).printAll("");
				
				for (int i = 0; i < xBins; i++){
					for (int j = 0; j < yBins; j++){
						for (int i2 = 0; i2 < xBins; i2++){
							for (int j2 = 0; j2 < yBins; j2++){
								int here1 = i*yBins + j;
								int here2 = i2*yBins + j2;
								int here = here1*xBins*yBins + here2;
								//Gamma[here] += lambdaG * G.xEigen[xEig][i+1] * G.yEigen[yEig][j+1] * G.xEigen[xEig][i2+1] * G.yEigen[yEig][j2+1];
								Gamma[here] += lambdaG * EigenVectorWithoutUnobservedBaslines.real[i*yBins + j] * EigenVectorWithoutUnobservedBaslines.real[i2*yBins + j2];
							}
						}
					}
				}
			}
		}
	}
	

	//Compute Eigensystem
	gsl_matrix_view m = gsl_matrix_view_array (Gamma, xBins*yBins, xBins*yBins);
	gsl_vector *eval = gsl_vector_alloc (xBins*yBins);
	gsl_matrix *evec = gsl_matrix_alloc (xBins*yBins, xBins*yBins);
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (xBins*yBins);
	gsl_eigen_symmv (&m.matrix, eval, evec, w);
	gsl_eigen_symmv_free (w);
	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
	
	int nEigToModify = 0;
	//Compute and store eigenvalues and fourier transformed eigenvectors	
	for (int n = 0; n < xBins*yBins; n++){
		double lambdaGamma = gsl_vector_get(eval, n);
		if (lambdaGamma / precon.noiseMin > PreconEVThreshold) {
			CVector eigenvector(s);
			int counter = 0;
			for (int i = 0; i < xBins; i++){
				for (int j = 0; j < yBins; j++){
					for (int k2 = 0; k2 < fBins; k2++){
						eigenvector.real[counter] = gsl_matrix_get(evec,i*yBins+j,n) * G.fdEigenMasked[k][k2+1];
						counter++;
					}
				}	
			}
			applyMask(eigenvector, precon.mask);
			CVector eigFourier = eigenvector.ijk2uvk();
			CVector eigFourierTrans = eigenvector.uvk2ijk();
			double barFactor = 0;
			for (int i = 0; i < nElements; i++){
				if (N.Cov.entry[i] != -1){
					barFactor += (eigFourier.real[i]*eigFourierTrans.real[i] - eigFourier.imag[i]*eigFourierTrans.imag[i]) / N.avgCov.entry[i];
				}
			}
			if (lambdaGamma * barFactor > PreconEVThreshold){
				precon.GammaEigenvectorsFourier[k].push_back(eigFourier);
				precon.GammaEigenvectorsFourierTrans[k].push_back(eigFourierTrans);
				precon.GammaEigenvalues[k].push_back(lambdaGamma);
				precon.GammaBarFactors[k].push_back(barFactor);
				nEigToModify++;
			}
		}
	}
	
	gsl_vector_free (eval);
	gsl_matrix_free (evec);	
	delete[] Gamma;
	return nEigToModify;
}

void constructPreconditioner(Preconditioner& precon, ResolvedCovariance& R, UnresolvedCovariance& U, GalacticCovariance& G, NoiseCovariance& N, CVector& mask){
	cout << "Now constructing the preconditioner..." << endl << endl;
	precon.mask = mask;
	
	//Calculate the mininum value of the noise, which gives the largest factor by which to multiply the eigenvalues in the new basis
	precon.noiseMin = 1e50;
	for (int n = 0; n < N.Cov.nElements; n++){
		if (N.Cov.entry[n] < precon.noiseMin && N.Cov.entry[n] != -1) precon.noiseMin = N.Cov.entry[n];
	}
	
	//Calculate the bar factor for R, which is always the same.
	double RbarFactor = 0;
	if (!removeR){
	 	cout << "Now determining which R eigenvectors to include..." << endl << endl;
		CVector sampleEVforRbar = CVector(s);
		for (int k = 0; k < fBins; k++) sampleEVforRbar.real[k] = R.EigenMasked[0][0][k+1];
		CVector FFTsampleEVforRbar = sampleEVforRbar.ijk2uvk();
		CVector IFFTsampleEVforRbar = sampleEVforRbar.uvk2ijk();
		for (int n = 0; n < sampleEVforRbar.nElements; n++){
			RbarFactor += (FFTsampleEVforRbar.real[n]*IFFTsampleEVforRbar.real[n] - FFTsampleEVforRbar.imag[n]*IFFTsampleEVforRbar.imag[n]) / N.avgCov.entry[n];
		} 
	}
	

	vector<double> biggestLambdaBar(G.fdEigen.size(), 0);
	vector<int> biggestLambdaBarIndex(G.fdEigen.size(), -1);
	vector<int> numberOfImportantEigenvalues(G.fdEigen.size(), 0);	
	
	int ReigsToModify = 0;
	int GeigsToModify = 0;
	int GammaEigsToModify = 0;
	
	//Figure out which eigenvectors of Rxy will be important.  Store them in R.eigenToUseInGamma[f Eigenvector number][which point sources]
	if (!removeR){
		for (int k = 0; k < G.fdEigen.size(); k++){
			vector<int> RtoUse;
			for (int n = 0; n < R.nSources; n++){
				double lambdaRBar = 0;
				if (R.EigenMasked[n].size() > k){
					lambdaRBar = RbarFactor * R.EigenMasked[n][k][0];
				}
				
				if (lambdaRBar > PreconEVThreshold) {
					RtoUse.push_back(n);
					numberOfImportantEigenvalues[k]++;
					ReigsToModify++;
				}
				if (lambdaRBar > biggestLambdaBar[k]){
					biggestLambdaBar[k] = lambdaRBar;
					biggestLambdaBarIndex[k] = n;
				}
			}
			if (RtoUse.empty()) RtoUse.push_back(-1); //A flag to indicate that no eigenvectors are important at this level.
			R.eigenToUseInGamma.push_back(RtoUse);					
		} 
	}
	
	//Calculate which eigenvectors of Gxy will be important.  Store them in G.eigenToUseInGamma[f Eigenvector][list of x-y pairs][0=x,1=y]
	cout << "Now determining which G eigenvectors to include..." << endl << endl;
	vector< vector< vector<double> > > GbarFactors(G.xEigen.size(), vector< vector<double> >(G.yEigen.size(), vector<double>(fBins,0)));
	for (int k = 0; k < G.fdEigenMasked.size(); k++){
		vector< vector<int> > eigenvectorsToUse;
		if (k >= G.fdEigenMasked.size()) {
			eigenvectorsToUse.push_back(vector<int>(2,-1));
		} else {
			for (int l = 0; l < G.xEigen.size(); l++){
				for (int m = 0; m < G.yEigen.size(); m++){
					double lambdaG = G.xEigen[l][0]*G.yEigen[m][0]*G.fdEigenMasked[k][0];
					if ((lambdaG / precon.noiseMin) > PreconEVThreshold){
						CVector v(s);
						int counter = 0;
						for (int i = 0; i < xBins; i++){
							for (int j = 0; j < yBins; j++){
								for (int k2 = 0; k2 < fBins; k2++){
									v.real[counter] = G.xEigen[l][i+1] * G.yEigen[m][j+1] * G.fdEigenMasked[k][k2+1];
									counter++;
								}
							}	
						}
						CVector vFourier = v.ijk2uvk();
						CVector vFourierTrans = v.uvk2ijk();
						double barFactor = 0;
						for (int i = 0; i < nElements; i++){
							if (N.Cov.entry[i] != -1){
								barFactor += (vFourier.real[i]*vFourierTrans.real[i] - vFourier.imag[i]*vFourierTrans.imag[i]) / N.avgCov.entry[i]; 
							}
						}
						double lambdaGbar = lambdaG*barFactor;
						GbarFactors[l][m][k] = barFactor;
						if (lambdaG*barFactor > PreconEVThreshold) {
							vector<int> coord(2,0);
							coord[0] = l;
							coord[1] = m;
							eigenvectorsToUse.push_back(coord);
							numberOfImportantEigenvalues[k]++;
							GeigsToModify++;
						}
						if (lambdaGbar > biggestLambdaBar[k]){
							biggestLambdaBar[k] = lambdaGbar;
							biggestLambdaBarIndex[k] = -1; //flag to use Geig with x=0, y=0
						}
					}
				}
			}
		}
		if (eigenvectorsToUse.empty()) eigenvectorsToUse.push_back(vector<int>(2,-1)); //flag to indicate that there are no eigenvectors to use for this f eigenvalue
		G.eigenToUseInGamma.push_back(eigenvectorsToUse);
	}
	
	if (biggestLambdaBar[0] > 1.0e13) {
		cout << endl <<  "************************************************************" << endl;
		cout << "WARNING: INSUFFIENT PRECISION!!!" << endl;
		cout << "************************************************************" << endl << endl;
	}
	
	//Find and store the eigenvalues and fourier transfored eigenvectors of Gamma 
	cout << "Now calculating the Gamma eigensystems..." << endl << endl;
	precon.GammaEigenvectorsFourier = vector< vector<CVector> >(G.fdEigenMasked.size(),vector<CVector>(0,CVector()));
	precon.GammaEigenvectorsFourierTrans = vector< vector<CVector> >(G.fdEigenMasked.size(),vector<CVector>(0,CVector()));
	precon.GammaEigenvalues = vector< vector<double> >(G.fdEigenMasked.size(),vector<double>(0,0));
	precon.GammaBarFactors = vector< vector<double> >(G.fdEigenMasked.size(),vector<double>(0,0));
	for (int k = 0; k < G.fdEigenMasked.size(); k++){	
		int nEigToModify = 0;
		if (((!removeR && R.eigenToUseInGamma[k][0] != -1) ) || (G.eigenToUseInGamma[k][0][0] != -1)){
			nEigToModify = GammaEigenDecomp(k,precon,R,G,N);
			GammaEigsToModify += nEigToModify;
		}
		cout << k+1 << " of " <<  G.fdEigenMasked.size() << " eigensystems done. " << nEigToModify << " eigenvalues of this Gamma need preconditioning." << endl;
	}
	cout << endl;
	
	int UeigsToModify = 0;
	for (int n = 0; n < U.fdEigenMasked.size(); n++){
		if (U.fdEigenMasked[n][0]/precon.noiseMin > PreconEVThreshold) UeigsToModify++;
	}
/* 	cout << nElements << endl;
	cout << ReigsToModify << endl;
	cout << UeigsToModify << endl;
	cout << GeigsToModify << endl;
	cout << GammaEigsToModify << endl; */
}

//Step 1: Apply N based preconditioner
void applyPN(CVector& Px, Preconditioner& precon, NoiseCovariance& N, bool transpose, CVector& wedgeMask){
	CVector PxFourier = Px;
	if (projectOutHorizonWedge && (CequalsN || CequalsI)){
		if (transpose){
			for (int n = 0; n < Px.nElements; n++){
				if (wedgeMask.real[n]){
					PxFourier.real[n] =  Px.real[n]/sqrt(N.avgCov.entry[n]);
					PxFourier.imag[n] =  Px.imag[n]/sqrt(N.avgCov.entry[n]);
				}
			}
			Px = PxFourier.uvw2ijk();
		} else {
			PxFourier = Px.ijk2uvw();
			for (int n = 0; n < Px.nElements; n++){
				if (wedgeMask.real[n]){
					Px.real[n] =  PxFourier.real[n]/sqrt(N.avgCov.entry[n]);
					Px.imag[n] =  PxFourier.imag[n]/sqrt(N.avgCov.entry[n]);
				}
			}
		}
	} else {
		if (transpose){
			for (int n = 0; n < Px.nElements; n++){
				if (CequalsN || CequalsI){
					PxFourier.real[n] =  Px.real[n]/sqrt(N.effCov.entry[n]);
					PxFourier.imag[n] =  Px.imag[n]/sqrt(N.effCov.entry[n]);
				} else {
					PxFourier.real[n] =  Px.real[n]/sqrt(N.avgCov.entry[n]);
					PxFourier.imag[n] =  Px.imag[n]/sqrt(N.avgCov.entry[n]);
				}
			}
			Px = PxFourier.uvk2ijk();
		} else {
			PxFourier = Px.ijk2uvk();
			for (int n = 0; n < Px.nElements; n++){
				if (CequalsN || CequalsI){
					Px.real[n] =  PxFourier.real[n]/sqrt(N.effCov.entry[n]);
					Px.imag[n] =  PxFourier.imag[n]/sqrt(N.effCov.entry[n]);  //TODO: Should this be eff cov??????
				} else {
					Px.real[n] =  PxFourier.real[n]/sqrt(N.avgCov.entry[n]);
					Px.imag[n] =  PxFourier.imag[n]/sqrt(N.avgCov.entry[n]); 
				}
			}
		} //Note: for the modes where we don't want to act for the pseudoinverse, this simply divides by 1 (based on how effCov is defined) 
	}
}

//Step 2: Apply R+G based preconditioner
void applyPGamma(CVector& Px, Preconditioner& precon, UnresolvedCovariance& U, GalacticCovariance& G, NoiseCovariance& N, bool transpose){
	CVector PxCopy = Px;
	CVector NFourierTimesPx = Px;
	//Multiply by Nfourier to the +/- 1
	if (transpose){
		for (int n = 0; n < nElements; n++) NFourierTimesPx.real[n] /= sqrt(N.avgCov.entry[n]);
		for (int n = 0; n < nElements; n++) NFourierTimesPx.imag[n] /= sqrt(N.avgCov.entry[n]);
	} else {
		for (int n = 0; n < nElements; n++) NFourierTimesPx.real[n] *= sqrt(N.avgCov.entry[n]);
		for (int n = 0; n < nElements; n++) NFourierTimesPx.imag[n] *= sqrt(N.avgCov.entry[n]);
	} //Note: for the modes where we don't want to act for the pseudoinverse, this simply divides by 1 (based on how effCov is defined) 
	
	
	for (int k = 0; k < G.fdEigen.size(); k++){
		if (precon.GammaEigenvectorsFourier.size() > k){
			if (!(precon.GammaEigenvectorsFourier[k].empty())){
				for (int n = 0; n < precon.GammaEigenvectorsFourier[k].size(); n++){
					double lambdaUbar =  U.fdEigenMasked[k][0] * precon.GammaBarFactors[k][n];
					double lambdaGammabar = precon.GammaEigenvalues[k][n] * precon.GammaBarFactors[k][n];
					double beta = 1 - (1 + lambdaUbar)/sqrt((1 + lambdaUbar)*(1 + lambdaUbar + lambdaGammabar));
					
					//Compute Px dotted into vFourierTranspose
					double vDotPxReal = 0;
					double vDotPxImag = 0;
					for (int i = 0; i < xBins*yBins; i++){
						for (int j = 0; j < fBins; j++){
							int here = i*fBins + j;
							if (N.Cov.entry[here] != -1){
								vDotPxReal += precon.GammaEigenvectorsFourierTrans[k][n].real[here] * NFourierTimesPx.real[here] - precon.GammaEigenvectorsFourierTrans[k][n].imag[here] * NFourierTimesPx.imag[here];
								vDotPxImag += precon.GammaEigenvectorsFourierTrans[k][n].imag[here] * NFourierTimesPx.real[here] + precon.GammaEigenvectorsFourierTrans[k][n].real[here] * NFourierTimesPx.imag[here];
							}
						}
					}
					
					//Multiply that dot product by vFourier
					for (int i = 0; i < xBins*yBins; i++){
						for (int j = 0; j < fBins; j++){
							int here = i*fBins + j;
							if (N.Cov.entry[here] != -1){
								double productReal = vDotPxReal*precon.GammaEigenvectorsFourier[k][n].real[here] + vDotPxImag*precon.GammaEigenvectorsFourier[k][n].imag[here];
								double productImag = vDotPxImag*precon.GammaEigenvectorsFourier[k][n].real[here] + vDotPxReal*precon.GammaEigenvectorsFourier[k][n].imag[here];
								if (transpose){
									PxCopy.real[here] -= beta * sqrt(N.avgCov.entry[here]) * productReal;
									PxCopy.imag[here] -= beta * sqrt(N.avgCov.entry[here]) * productImag;
								} else {
									PxCopy.real[here] -= beta / sqrt(N.avgCov.entry[here]) * productReal;
									PxCopy.imag[here] -= beta / sqrt(N.avgCov.entry[here]) * productImag;
								}
							}
						}
					}
				}
			}
		}
	}
						
	Px = PxCopy;	
}

//Step 3: Apply U based preconditioner
void applyPU(CVector& Px, Preconditioner& precon, UnresolvedCovariance& U, NoiseCovariance& N){	
	CVector PxCopy = Px;
	for (int n = 0; n < U.fdEigenMasked.size(); n++){
		double lambdaU = U.fdEigenMasked[n][0];
		if (lambdaU/precon.noiseMin > PreconEVThreshold){
			vector<double> v(fBins,0);
			for (int k = 0; k < fBins; k++) v[k] = U.fdEigenMasked[n][k+1];
			for (int i = 0; i < xBins; i++){
				for (int j = 0; j < yBins; j++){
					int thisBlock = i*fBins*yBins + j*fBins;
					if (N.Cov.entry[thisBlock] != -1){
						double beta = 1 - sqrt(1/(1+lambdaU/N.avgCov.entry[thisBlock]));
						double vDotPxReal = 0;
						double vDotPxImag = 0;
						for (int k = 0; k < fBins; k++){
							vDotPxReal += v[k]*Px.real[thisBlock + k];
							vDotPxImag += v[k]*Px.imag[thisBlock + k];
						}
						for (int k = 0; k < fBins; k++){
							PxCopy.real[thisBlock + k] -= beta*(v[k]*vDotPxReal);
							PxCopy.imag[thisBlock + k] -= beta*(v[k]*vDotPxImag);
						}
					}
				}
			}
		}
	}
	Px = PxCopy;
}

//Wrapper function for applying all the different preconditioners
CVector precondition(CVector& x, Preconditioner& precon, ResolvedCovariance& R, UnresolvedCovariance& U, GalacticCovariance& G, NoiseCovariance& N, bool transpose, CVector& wedgeMask){
	if (turnOffPreconditioner) return x;
	CVector Px = x;
	
	if (CequalsN || CequalsI){
		if (projectOutHorizonWedge){
			if (transpose) {
				Px = Px.ijk2uvw();
				applyPN(Px,precon,N,transpose,wedgeMask);
			} else {
				applyPN(Px,precon,N,transpose,wedgeMask);
				Px = Px.uvw2ijk();
			}
		} else {
			if (transpose) {
				Px = Px.ijk2uvk();
				applyPN(Px,precon,N,transpose,wedgeMask);
			} else {
				applyPN(Px,precon,N,transpose,wedgeMask);
				Px = Px.uvk2ijk();
			}
		}
		for (int n = 0; n < nElements; n++) Px.imag[n] = 0;
		return Px;
	}
	
	if (CequalsNU){
		if (transpose) {
			Px = Px.ijk2uvk();
			applyPU(Px,precon,U,N);
			applyPN(Px,precon,N,transpose,wedgeMask);
		} else {
			applyPN(Px,precon,N,transpose,wedgeMask);
			applyPU(Px,precon,U,N);
			Px = Px.uvk2ijk();
		}
		return Px;
	}
	
	if (transpose) {
		Px = Px.ijk2uvk();
		applyPU(Px,precon,U,N);
		applyPGamma(Px,precon,U,G,N,transpose);
		applyPN(Px,precon,N,transpose,wedgeMask);
	} else {
		applyPN(Px,precon,N,transpose,wedgeMask);
		applyPGamma(Px,precon,U,G,N,transpose);
		applyPU(Px,precon,U,N);
		Px = Px.uvk2ijk();
	}

	for (int n = 0; n < nElements; n++) Px.imag[n] = 0;
	return Px;
	
}

/*******************************************************************
GENERATING THE SIMULATED DATA CUBE
********************************************************************/

//Given a list of point sources with given properties, generates a contribution to x from resolved point sources that reasonably agrees with the believed value
CVector generatexR(ResolvedCovariance& R) {
	CVector xR = CVector(s);  //Technique using random fluxes and spectral indices.
	for (int n = 0; n < R.nSources; n++){ //Old Spectral Index Method
		double r1 = 1.0*rand()/RAND_MAX;
		double r2 = 1.0*rand()/RAND_MAX;
		double flux = R.fluxes[n]*(1 + R.fluxPercentError*sqrt(-2.0*log(r1))*cos(2.0*pi*r2));
		double spectralIndex = R.spectralIndices[n]*(1 + R.spectralIndexPercentError*sqrt(-2.0*log(r1))*sin(2.0*pi*r2));
		for (int k = 0; k < fBins; k++){
			xR.real[R.positions[n][0]*xR.fBins*xR.yBins + R.positions[n][1]*xR.fBins + k] = 1.4e-3 * flux / R.omegaPix * pow(R.freqs[k]/R.referenceFrequency,-2-spectralIndex);
		}
	}
	if (removeMean) for (int n = 0; n < R.Dmean.nElements; n++) xR.real[n] -= R.Dmean.entry[n];
	
	return xR;
}

//Generates a contribution to x from unresolved point sources using prewhitening and Toeplitz matrix techniques
CVector generatexU(UnresolvedCovariance& U){

	vector< vector<double> > fluxRandomField = U.xCov.gaussRandField2D(U.yCov);  
	CVector xU(s);
	
	for (int n = 0; n < U.fEigen.size(); n++){
		double r1 = 1.0*rand()/RAND_MAX;
		double r2 = 1.0*rand()/RAND_MAX;
		for (int m = 0; m < fBins; m++){
			for (int i = 0; i < xBins; i++){
				for (int j = 0; j < yBins; j++){
					double temp = U.Dcross.entry[i*fBins*yBins + j*fBins + m] * fluxRandomField[i][j] * sqrt(-2.0*log(r1))*cos(2.0*pi*r2) * sqrt(U.fEigen[n][0]) * U.fEigen[n][m+1];
					xU.real[i*yBins*fBins + j*fBins + m] += temp;
				}
			}
		}
	}
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			for (int k = 0; k < fBins; k++){
				if (!removeMean) xU.real[i*yBins*fBins + j*fBins + k] += U.Dmean.entry[i*yBins*fBins + j*fBins + k];
			}
		}
	}
	
	return xU;
}

//Generates a contribution to x from the galactic synchrotron using prewhitening and Toeplitz matrix techniques
CVector generatexG(GalacticCovariance& G){

	vector< vector<double> > randomField(xBins, vector<double> (yBins,0));
	for (int n = 0; n < G.xEigen.size(); n++){
		for (int m = 0; m < G.yEigen.size(); m++){
			double r1 = 1.0*rand()/RAND_MAX;
			double r2 = 1.0*rand()/RAND_MAX;
			double gRand = sqrt(G.xEigen[n][0] * G.yEigen[m][0]) * sqrt(-2.0*log(r1))*cos(2.0*pi*r2);
			for (int i = 0; i < xBins; i++){
				for (int j = 0; j < yBins; j++){
					randomField[i][j] += gRand * G.xEigen[n][i+1] * G.yEigen[m][j+1];
				}
			}
		}
	}

	CVector xG(s);
	for (int n = 0; n < G.fEigen.size(); n++){
		double r1 = 1.0*rand()/RAND_MAX;
		double r2 = 1.0*rand()/RAND_MAX;
		for (int m = 0; m < fBins; m++){
			for (int i = 0; i < xBins; i++){
				for (int j = 0; j < yBins; j++){
					double temp = G.Dcross.entry[i*fBins*yBins + j*fBins + m] * randomField[i][j] * sqrt(-2.0*log(r1))*cos(2.0*pi*r2) * sqrt(G.fEigen[n][0]) * G.fEigen[n][m+1];
					xG.real[i*yBins*fBins + j*fBins + m] += temp;
				}
			}
		}
	}
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			for (int k = 0; k < fBins; k++){
				if (!removeMean) xG.real[i*yBins*fBins + j*fBins + k] += G.Dmean.entry[i*yBins*fBins + j*fBins + k];
			}
		}
	}
	
	return xG;
}


//Given an N covariance matrix in the uvk basis, this function creates a zero-mean signal with this covariance
CVector generatexN(NoiseCovariance& N){ 
	CVector xNFourier = CVector(s);
	RDMatrix used = RDMatrix(s);
	for (int n = 0; n < xBins; n++){
		for (int m = 0; m < yBins; m++){
			for (int l = 0; l < fBins; l++){
				int here = l + m*fBins + n*fBins*yBins;
				if (used.entry[here] == 0){ //if this place wasn't already the complex conjugate of somewhere else
					used.entry[here] = 1;
					double r1 = 1.0*rand()/RAND_MAX;
					double r2 = 1.0*rand()/RAND_MAX;
					double gRandReal = 0;
					double gRandImag = 0;
					if (N.Cov.entry[here] != -1) gRandReal = sqrt(N.effCov.entry[here])*sqrt(-2.0*log(r1))*cos(2.0*pi*r2); //gaussian random variable
					xNFourier.real[here] = gRandReal;
					if (!((n == 0 || n == xBins/2) && (m == 0 || m == yBins/2) /*&& (l == 0 || l == fBins/2)*/)){ //If the entry here needs to be complex
						if (N.Cov.entry[here] != -1) gRandImag = sqrt(N.effCov.entry[here])*sqrt(-2.0*log(r1))*sin(2.0*pi*r2);
						xNFourier.real[here] = gRandReal/pow(2,.5);
						xNFourier.imag[here] = gRandImag/pow(2,.5);
						int lOpposite = l; //calculates the location of the opposite
						int mOpposite = xNFourier.yBins - m;
						int nOpposite = xNFourier.xBins - n;
						if (l == 0) lOpposite = l;
						if (m == 0) mOpposite = m;
						if (n == 0) nOpposite = n;
						int opposite = lOpposite + mOpposite*xNFourier.fBins + nOpposite*xNFourier.fBins*xNFourier.yBins;
						if (used.entry[opposite] == 0){
							if (fabs(N.effCov.entry[here] - N.effCov.entry[opposite]) > 1e-12*N.effCov.entry[opposite]) cout << "(" << n << "," << m << "," << l << ") to (" << nOpposite << "," << mOpposite << "," << lOpposite << "): N here and opposite: " << N.effCov.entry[here] << "  " << N.effCov.entry[opposite] << endl;							
							used.entry[opposite] = 1;
							if (N.effCov.entry[here] != -1){
								xNFourier.real[opposite] = gRandReal/pow(2,.5); //complex conjugate
								xNFourier.imag[opposite] = -gRandImag/pow(2,.5);
							}
						}
					} 
				}
			}
		}
	}
	
	CVector xN = xNFourier.uvk2ijk();
	bool warning = false;
	for (int n = 0; n < xN.nElements; n++){
		if (!warning && fabs(xN.imag[n]) > 1e-9){
			cout << endl << "WARNING: imaginary value of " << xN.imag[n] << " detected in xN at n = " << n << " (and set to 0)!" << endl << endl;
			warning = true;
		}
		xN.imag[n] = 0;
		xN.real[n] *= sqrt(xBins * yBins);
	}

	double rms = 0;
	for (int n = 0; n < xN.nElements; n++) rms += pow(xN.real[n],2)/xN.nElements;
	rms = sqrt(rms);
	N.Cov.rms = rms;
	if(verbose) cout << "Noise RMS = " << rms <<  " K" << endl;
	return xN;
}

//Returns one instantiation of x, making sure that it is zeroPadded (all zeros outside the region of interest)
CVector generateUniverse(ResolvedCovariance& R, UnresolvedCovariance& U, GalacticCovariance& G, NoiseCovariance& N){
	if (verbose) cout << "Now generating an artificial data cube..." << endl << endl;	
	if (CequalsN || CequalsI) return generatexN(N);	
	
	if (CequalsNU){
		CVector xU = generatexU(U);
		CVector xN = generatexN(N);
		return xU + xN;
	}
	
	CVector xR(s);
	CVector xG(s);
	if (!removeR) xR = generatexR(R);
	CVector xU = generatexU(U);
	if (!removeG) xG = generatexG(G);
	CVector xN = generatexN(N);
	if (printDataCubesToFile){
		xR.printRealToFile("xR.dat");
		xU.printRealToFile("xU.dat");
		xG.printRealToFile("xG.dat");
		xN.printRealToFile("xN.dat");
	}
	return xR + xU + xG + xN; 
}

//Testing to see if the covariance of the random universes behaves properly
void CovarianceTesting(ResolvedCovariance& R, UnresolvedCovariance& U, GalacticCovariance& G, NoiseCovariance& N, string matrixName){
	vector< vector<double> > second(xBins*yBins*fBins, vector<double>(xBins*yBins*fBins,0));
	vector<double> average(xBins*yBins*fBins,0);
	
	int trials = covTestingTrials;
	cout << endl << "Now testing the covariance of Generate" << matrixName << "..." << endl;
	int percent = 0;
	
	for (int n = 0; n < trials; n++){
		if (100.0*n/trials > percent) {
			percent = int(ceil(100.0*n/trials));
			cout << percent << "% complete." << endl;
		}
		CVector x(s);
		if (matrixName.compare("R") == 0) x = generatexR(R);
		else if (matrixName.compare("U") == 0) x = generatexU(U);
		else if (matrixName.compare("G") == 0) x = generatexG(G);
		else if (matrixName.compare("N") == 0) x = generatexN(N);
		else if (matrixName.compare("C") == 0) x = generateUniverse(R,U,G,N);
		else {
			cout << endl << "WARNING: Trying test the covariance of a matrix that doesn't exist." << endl << endl;
			return;
		}
		for (int j = 0; j < xBins*yBins*fBins; j++){ 
		average[j] += x.real[j]/trials;
			for (int k = 0; k < xBins*yBins*fBins; k++){
				second[j][k] += x.real[j]*x.real[k]/trials;
			}
		}
	}
	cout << "Done.                  " << endl;
	
	ofstream outfile;
	string outFilename;
	if (matrixName.compare("R") == 0) outFilename = "RcovTest.dat";
	if (matrixName.compare("U") == 0) outFilename = "UcovTest.dat";
	if (matrixName.compare("G") == 0) outFilename = "GcovTest.dat";
	if (matrixName.compare("N") == 0) outFilename = "NcovTest.dat";
	if (matrixName.compare("C") == 0) outFilename = "CcovTest.dat";
	outfile.open(outFilename.c_str(), ios::trunc);	
	outfile.precision(30);
	for (int i = 0; i < xBins*yBins*fBins; i++){
		for (int j = 0; j < xBins*yBins*fBins; j++){
			outfile << (second[i][j] - average[i]*average[j])*(trials/(trials-1)) << ", ";
		}
		outfile << endl;
	}
	outfile.close();
	
	return;
}


/*******************************************************************
MULTIPLY BY THE COVARIANCE MATRICES
********************************************************************/

CVector multiplyByR(CVector& y, ResolvedCovariance& R){
	CVector yPrewhiteneDmean = R.Dmean*y;
	CVector RMeany(s);
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			double sum = 0;
			for (int k = 0; k < fBins; k++) sum += yPrewhiteneDmean.real[i*fBins*yBins + j*fBins + k];
			for (int k = 0; k < fBins; k++) RMeany.real[i*fBins*yBins + j*fBins + k] = sum;
		}
	}
	RMeany = R.Dmean * RMeany;
	
	CVector RCrossy(s);
	RCrossy = R.Dcross * y;
	for (int n = 0; n < R.nSources; n++){
		RCrossy = R.Cov[n] * RCrossy;		
	}
	RCrossy = R.Dcross * RCrossy;
	
	return RCrossy - RMeany;
}

CVector multiplyByU(CVector& y, UnresolvedCovariance& U){
	return U.Dcross * (U.xCov * (U.yCov * (U.fCov * (U.Dcross * y))));
}

CVector multiplyByG(CVector& y,GalacticCovariance& G){
	return G.Dcross * (G.xCov * (G.yCov * (G.fCov * (G.Dcross * y))));
}

CVector multiplyByN(CVector& y, NoiseCovariance& N){
	CVector Ny = (N.effCov*(y.ijk2uvk())).uvk2ijk(); //Note that I use the effective covariance here, which has all the flags replaced with the largest number
	bool warned = false;
	for (int n = 0; n < Ny.nElements; n++){
		if (!warned && fabs(Ny.imag[n]) > 1e-10){
			cout << endl << "WARNING: imaginary value of " << Ny.imag[n] << " multiplication by N at n = " << n << "(and set to zero)." << endl << endl;
			warned = true;
		}
		Ny.imag[n] = 0;
	}	
	return Ny;
}


CVector multiplyByC(CVector& x, ResolvedCovariance& R, UnresolvedCovariance& U, GalacticCovariance& G, NoiseCovariance& N){
	if (CequalsN || CequalsI) return multiplyByN(x,N);
	
	if (CequalsNU){
		CVector Cx(s);
		Cx = Cx + multiplyByU(x,U);
		Cx = Cx + multiplyByN(x,N);
		return Cx;
	}
	
	CVector Cx(s);
	if (!removeR) Cx = Cx + multiplyByR(x,R);
	Cx = Cx + multiplyByU(x,U);
	if (!removeG) Cx = Cx + multiplyByG(x,G);
	Cx = Cx + multiplyByN(x,N);
	return Cx;
}

void zeroOutUnobservedBaselines(CVector& dataCube, NoiseCovariance& N){
	CVector xFourier = dataCube.ijk2uvk();
	for (int n = 0; n < xBins*yBins*fBins; n++){ //Downweight all power in unobserved baselines
		if (N.Cov.entry[n] == -1){
			xFourier.real[n] = 0;
			xFourier.imag[n] = 0;
		}
	}
	dataCube = xFourier.uvk2ijk();
}

CVector calculateWedgeMask(){
	double comovingDist = 0;
	double zEnd = f21cm/fStart - 1;
	double deltaRedshift = .00001;
	double zLeft = 0;
	double zRight = deltaRedshift;
	while (true){
		if ((zLeft + zRight)/2.0 > zEnd) break;
		comovingDist += c/H0*((1.0/sqrt(OmegaM*pow(1+zRight,3) + OmegaL) + 4.0/sqrt(OmegaM*pow(1+(zLeft + zRight)/2,3) + OmegaL) + 1.0/sqrt(OmegaM*pow(1+zLeft,3)+OmegaL))*deltaRedshift/6);		
		zLeft = zRight;
		zRight += deltaRedshift;
	}
	comovingDist -= fLength/2.0;
	double zCenter = (zLeft + zRight) / 2;
	double wedgeCoeff = sqrt(OmegaM*pow((1+zCenter),3)+OmegaL)*(comovingDist)/(c/H0)/(1+zCenter);

	CVector wedgeMask(s);
	double deltaKPara = 2.0*pi / fLength;
	double deltaKPerp = 2.0*pi / xyLength;
	for (int k = 0; k < fBins; k++){
		double kParaHere = fabs(k-fBins/2)*deltaKPara;
		for (int i = 0; i < xBins; i++){
			for (int j = 0; j < yBins; j++){
				double kPerpHere = sqrt(pow((i-xBins/2),2) + pow((j-yBins/2),2))*deltaKPerp;
				if (kParaHere > (wedgeCoeff * kPerpHere + horizonWedgeBuffer)){
					wedgeMask.real[i*yBins*fBins + j*fBins + k] = 1;
				} else {
					wedgeMask.real[i*yBins*fBins + j*fBins + k] = 0;
				}
			}
		}
	}
	return wedgeMask;
}

void projectOutHorizonWedgePlusBuffer(CVector& x, CVector& wedgeMask){
	CVector xFT = x.ijk2uvw();
	CVector xFiltered = xFT;
	for (int n = 0; n < nElements; n++){
		xFiltered.real[n] *= wedgeMask.real[n];
		xFiltered.imag[n] *= wedgeMask.real[n];
	}
	x = xFiltered.uvw2ijk();
}

void projectOutMaskedOrUnobservedModes(CVector& x, CVector& mask, NoiseCovariance& N, CVector& wedgeMask){
	applyMask(x,mask);
	zeroOutUnobservedBaselines(x,N);
	if (projectOutHorizonWedge) projectOutHorizonWedgePlusBuffer(x,wedgeMask);
}

CVector multiplyByPsuedoC(CVector& x, ResolvedCovariance& R, UnresolvedCovariance& U, GalacticCovariance& G, NoiseCovariance& N, CVector& mask, CVector& wedgeMask){
	CVector projectedX = x;
	projectOutMaskedOrUnobservedModes(projectedX, mask, N, wedgeMask);
	
	CVector PiCPix = multiplyByC(projectedX,R,U,G,N);
	projectOutMaskedOrUnobservedModes(PiCPix, mask, N, wedgeMask);
	
	return PiCPix + x - projectedX;
}


CVector multiplyByPCPtrans(CVector& x, ResolvedCovariance& R, UnresolvedCovariance& U, GalacticCovariance& G, NoiseCovariance& N, Preconditioner& precon, CVector& mask, CVector& wedgeMask){
	CVector result = precondition(x,precon,R,U,G,N,true,wedgeMask);
	result = multiplyByPsuedoC(result,R,U,G,N,mask,wedgeMask);
	result = precondition(result,precon,R,U,G,N,false,wedgeMask);
	return result;	
}

void MatrixReconstruction(ResolvedCovariance& R, UnresolvedCovariance& U, GalacticCovariance& G, NoiseCovariance& N, CVector& mask){
	vector< vector<double> > reconstructed(xBins*yBins*fBins, vector<double>(xBins*yBins*fBins,0));
	int m = 0;
	cout << endl << "Now reconstructing C from basis vectors..." << endl;
	
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			cout << " " << floor(100.0 * i/xBins + 100.0*j/xBins/yBins) << "% done. \r" << std::flush;
			for (int k = 0; k < fBins; k++){
				CVector unitVector(s);
				unitVector.real[i*yBins*fBins + j*fBins + k] = 1;
				CVector CTimesUnit(s);
				CTimesUnit = multiplyByC(unitVector,R,U,G,N);				
				for (int i2 = 0; i2 < xBins; i2++){
					for (int j2 = 0; j2 < yBins; j2++){
						for (int k2 = 0; k2 < fBins; k2++){
							reconstructed[m][i2*yBins*fBins + j2*fBins + k2] = CTimesUnit.real[i2*yBins*fBins + j2*fBins + k2];
						}
					}
				}
				m++;
			}
		}
	}
	cout << "Done.                  " << endl;  
	
	CVector randTest(s);
	for (int n = 0; n < nElements; n++){
		randTest.real[n] = 1.0*rand()/RAND_MAX;
	}
	
	CVector product1 = multiplyByC(randTest,R,U,G,N);
	CVector product2(s);
	for (int n = 0; n < nElements; n++){
		for (int m = 0; m < nElements; m++){
			product2.real[n] += reconstructed[n][m] * randTest.real[m];
		}
	}
	
	//(product1 - product2).printAll("");
	
	
	
	string outFilename = "Crecon.dat";
	ofstream outfile;
	outfile.precision(30);
	outfile.open(outFilename.c_str(), ios::trunc);
 	/*for (int n = 0; n < xBins*yBins*fBins; n++){
		for (int m = 0; m < xBins*yBins*fBins; m++){
			outfile << reconstructed[m][n] << " ";
		}
		outfile << endl; 
	}*/
	
	//Output in a format for Adrian
	int nMasked = 0;
	for (int k = 0; k < fBins; k++) if (mask.real[k] == 1) nMasked++;
	//cut and reorder
	vector< vector<double> > cutAndReordered(xBins*yBins*(fBins-nMasked), vector<double>(xBins*yBins*(fBins-nMasked),0));
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			int fCounter1 = 0;
			for (int k = 0; k < fBins; k++){
				if (mask.real[i*yBins*fBins + j*fBins + k] == 0){
					for (int i2 = 0; i2 < xBins; i2++){
						for (int j2 = 0; j2 < yBins; j2++){
							int fCounter2 = 0;
							for (int k2 = 0; k2 < fBins; k2++){
								if (mask.real[i2*yBins*fBins + j2*fBins + k2] == 0){
									cutAndReordered[(fBins-nMasked - fCounter1 - 1)*xBins*yBins + j*xBins + (xBins-i-1)]
										[(fBins-nMasked - fCounter2 - 1)*xBins*yBins + j2*xBins + (xBins-i2-1)] =
										reconstructed[i*yBins*fBins + j*fBins + k][i2*yBins*fBins + j2*fBins + k2];
									fCounter2++;
								}
							}
						}
					}
					fCounter1++;
				}
			}
		}

	}
	
	for (int n = 0; n < (xBins*yBins*(fBins-nMasked)); n++){
		for (int m = 0; m <= n; m++){
			outfile << cutAndReordered[n][m] << endl;
		}
	}
	
	outfile.close();
}

/*******************************************************************
MAIN ALGORITHMS
********************************************************************/

CVector conjugateGradientAlgorithm(CVector& x, ResolvedCovariance& R, UnresolvedCovariance& U, GalacticCovariance& G, NoiseCovariance& N, Preconditioner& precon, CVector& mask, CVector& wedgeMask){
	if (verbose) cout << "Now starting Conjugate Gradient Algorithm..." << endl << endl;
	CVector xPrime = precondition(x,precon,R,U,G,N,false,wedgeMask);
	CVector yPrime_n(s);
	CVector y_np1(s);
	CVector rPrime_n = xPrime;
	CVector dPrime_n = rPrime_n;
	double convergenceStatistic_n = 1000000;
	for (int n = 0; n < x.nElements; n++){
		time_t start, end;
		time(&start);
		CVector PCPtdPrime_n = multiplyByPCPtrans(dPrime_n,R,U,G,N,precon,mask,wedgeMask);
		CVector alpha_n = (rPrime_n.dot(rPrime_n)) / (dPrime_n.dot(PCPtdPrime_n));
		//if (verbose) cout << "alpha_" << n << " = " << alpha_n.real[0] << " + " << alpha_n.imag[0] << "i." << endl;
		CVector yPrime_np1 = yPrime_n + (dPrime_n*alpha_n);
		y_np1 = precondition(yPrime_np1,precon,R,U,G,N,true,wedgeMask);
		for (int m = 0; m < y_np1.nElements; m++){
			if (fabs(y_np1.imag[m]) > 1e-10){
				cout << endl << "WARNING: imaginary value of " << y_np1.imag[m] << " detected in CGA at n = " << n << "!" << endl << endl;
				break;
			}
			y_np1.imag[m] = 0; //TODO: is this ok?
		}
		CVector convergenceTest = multiplyByPsuedoC(y_np1,R,U,G,N,mask,wedgeMask) - x;
		double convergenceStatistic_np1 = convergenceTest.magnitude()/x.magnitude();
		if (verbose) cout << n+1 << ") |(C*y - x)|/|x| = " << convergenceStatistic_np1 << endl;
		if (convergenceStatistic_np1 < CGBound){
			if (verbose) cout << endl << "The CGA has achieved an average magnitude of residual of less than " << CGBound << " in " << n+1 << " iteration(s)." << endl;
			if (verbose) cout << "Estimated Condition Number is less than " << pow( log(CGBound/2.0 * exp(2.0*(n+1))),2.0) << endl << endl;
			return y_np1;
		}
		if (fabs(convergenceStatistic_np1 - convergenceStatistic_n) < .00000000000001){
			cout << endl << "CGA has gotten stuck with convStat_old = " << convergenceStatistic_n << " and convStat_new = " << convergenceStatistic_np1 << "and is quitting." << endl;
			return y_np1;
		}
		convergenceStatistic_n = convergenceStatistic_np1;
		CVector rPrime_np1 = rPrime_n - (PCPtdPrime_n)*alpha_n;
		CVector beta_n = (rPrime_np1.dot(rPrime_np1)) / (rPrime_n.dot(rPrime_n));
		//if (verbose) cout << "beta_" << n << " = " << beta_n.real[0] << " + " << beta_n.imag[0] << "i." << endl;
		dPrime_n = rPrime_np1 + dPrime_n*beta_n;
		rPrime_n = rPrime_np1;
		yPrime_n = yPrime_np1;
		time(&end);
		cout << "       Calcuation Time: " << difftime (end,start) << " seconds." << endl;
		fftw_cleanup();
	}
	return y_np1;
}

void subtractMean(CVector& dataCube, ResolvedCovariance& R, UnresolvedCovariance& U, GalacticCovariance& G){
	if (CequalsN || CequalsI) return;

	if (removeMean) {
		for (int n = 0; n < dataCube.nElements; n++) dataCube.real[n] -= U.Dmean.entry[n];
		if (CequalsNU) return;
		if (!removeG){
			for (int n = 0; n < dataCube.nElements; n++) dataCube.real[n] -= G.Dmean.entry[n];	
		} 
		if (!removeR){
			for (int n = 0; n < dataCube.nElements; n++) dataCube.real[n] -= R.Dmean.entry[n];
		}
	}
}

CVector bandPowerSpectrum(CVector& yLeft, CVector& yRight, NoiseCovariance& N, bool useNoiseToFilter, vector<double>& sphericalPower){ //TODO: fix the fact that if xbins != ybins, the paraDist calculation is totally wrong.	
	int zxBins = zeroPad * xBins;
	int zyBins = zeroPad * yBins;
	int zfBins = zeroPad * fBins;
	
	bool useBothVectors = true;
	if (yRight.nElements == 1) useBothVectors = false;
				
	CVector yFourierSq(s,zxBins*zyBins*zfBins);
	if (useBothVectors){
		CVector yLeftZPFourier = yLeft.zeroPadAndFTforPowerSpectrum(false);
		CVector yRightZPFourier = yRight.zeroPadAndFTforPowerSpectrum(true);	
		for (int n = 0; n < yFourierSq.nElements; n++){ //The modulus squared of the 3D zero-padded Fourier transform
			yFourierSq.real[n] = (yLeftZPFourier.real[n]*yRightZPFourier.real[n] - yLeftZPFourier.imag[n]*yRightZPFourier.imag[n]); 
			yFourierSq.imag[n] = (yLeftZPFourier.real[n]*yRightZPFourier.imag[n] + yLeftZPFourier.imag[n]*yRightZPFourier.real[n]);
		}
		
	} else {
		CVector yZPFourier = yLeft.zeroPadAndFTforPowerSpectrum(false);
		for (int n = 0; n < yFourierSq.nElements; n++){ //The modulus squared of the 3D zero-padded Fourier transform
			yFourierSq.real[n] = (pow(yZPFourier.real[n],2) + pow(yZPFourier.imag[n],2)); 
			yFourierSq.imag[n] = 0;
		}
	}
	
	double dx = s->xyLength / xBins;
	double dy = s->xyLength / yBins;
	double df = s->fLength / fBins;
	double dkx = 2*pi / (zeroPad * s->xyLength);
	double dky = 2*pi / (zeroPad * s->xyLength);
	double dkf = 2*pi / (zeroPad * s->fLength);	
	vector<double> kX(zxBins,0);
	vector<double> kY(zyBins,0);
	vector<double> kF(zfBins,0);
	for (int i = -zxBins/2; i < zxBins/2; i++) kX[i + zxBins/2] = dkx*(i);
	for (int j = -zyBins/2; j < zyBins/2; j++) kY[j + zyBins/2] = dky*(j);
	for (int k = -zfBins/2; k < zfBins/2; k++) kF[k + zfBins/2] = dkf*(k);

	vector<int> perpCount(kPerpBins,0);
	//vector<vector<int> > perpTest(zxBins,vector<int>(zyBins,0));
	//cout << kParaBins << " " << kPerpBins << endl;
	CVector P(s,kParaBins*kPerpBins);
	for (int w = 0; w < zfBins; w++){
		int paraDist = int(round(abs(w-zfBins/2)/1.0/zeroPad));
		//if (paraDist < 0) paraDist = 0;
		double jf = 1;
		if (!JEquals1 && kF[w] != 0) jf = sin(kF[w]*df/2) / (kF[w]*df/2);
		for (int v = 0; v < zyBins; v++){
			double jy =  1;
			if (!JEquals1 && kY[v] != 0) jy = sin(kY[v]*dy/2) / (kY[v]*dy/2);
			for (int u = 0; u < zxBins; u++){ 
				int kSphereDist = int(round((sqrt(kF[w]*kF[w] + kY[v]*kY[v] + kX[u]*kX[u])/(2*pi / (s->fLength)))));
				int perpDist = int(round(sqrt(pow(u - zxBins/2,2) + pow(v - zyBins/2,2))/1.0/zeroPad));
				//if (perpDist < 0) perpDist = 0;
				double jx = 1;
				if (!JEquals1 && kX[u] != 0) jx = sin(kX[u]*dx/2) / (kX[u]*dx/2);	
				double normalization = pow(jx*jy*jf,2) / (s->xyLength * s->xyLength * s->fLength * zeroPad * zeroPad * zeroPad);
				double sphericalK = sqrt(pow(kX[u],2)+pow(kY[v],2)+pow(kF[w],2));
				if (!(w == zfBins/2 || (u == zxBins/2 && v == zyBins/2))){
					if (paraDist < kParaBins && perpDist < kPerpBins && paraDist >= 0 && perpDist >= 0){
						P.real[paraDist + perpDist*kParaBins] += .5 * normalization * yFourierSq.real[u*zyBins*zfBins + v*zfBins + w];
						perpCount[perpDist] += 1;
					}
					if (kSphereDist < kSphereBins && kSphereDist >= 0){
						sphericalPower[kSphereDist] += .5 * normalization * yFourierSq.real[u*zyBins*zfBins + v*zfBins + w];
					}
				}
			}
		}
	}

	return P;
}

void saveKBins(){
	//Save parallel 
	double deltaKPara = 2*pi / (s->fLength);
	double deltaKPerp = 2*pi / (s->xyLength); //TODO: This assumes that xBins == yBins and xLength == yLength
	double zeroBinWidthPara = pi/ (s->fLength) / (s->zeroPad);
	double zeroBinWidthPerp = pi/ (s->xyLength) / (s->zeroPad);
	ofstream outfile1, outfile2, outfile3;
	outfile1.open(kParaBinCentersFilename.c_str(), ios::trunc);
	outfile2.open(kPerpBinCentersFilename.c_str(), ios::trunc);
	
	outfile1 << (zeroBinWidthPara + deltaKPara*.5)*.5 << endl;
	for (int n = 0; n < kParaBins-1; n++){
		outfile1 << deltaKPara*(n+1) << endl;
	}
	outfile2 << 2.0/3.0 * (.5*deltaKPerp + pow(zeroBinWidthPerp,2) / (zeroBinWidthPerp + .5*deltaKPerp)) << endl;
	for (int n = 0; n < kPerpBins-1; n++){
		double innerRad = deltaKPerp*(n+.5);
		double outerRad = deltaKPerp*(n+1.5);		
		outfile2 << 2.0/3.0 * (outerRad + innerRad*innerRad / (innerRad + outerRad)) << endl;
	}
	outfile1.close();
	outfile2.close();
	double zeroSphericalBinWidth = sqrt(pow(pi/zeroPad/xyLength,2)+pow(pi/zeroPad/fLength,2));
	double deltaKSpherical =  2*pi / (s->fLength);
	outfile3.open(kSphericalBinCentersFilename.c_str(), ios::trunc);
	outfile3 << (zeroSphericalBinWidth + deltaKSpherical*.5)*.5 << endl;
	for (int n = 0; n < kParaBins-1; n++){
		outfile3 << deltaKSpherical*(n+1) << endl;
	}
	outfile3.close();

	if (saveKBinEdges){
		ofstream outfile4, outfile5, outfile6;
		outfile4.open(kParaBinEdgesFilename.c_str(), ios::trunc);
		outfile5.open(kPerpBinEdgesFilename.c_str(), ios::trunc);
		outfile6.open(kSphericalBinEdgesFilename.c_str(), ios::trunc);
		outfile4 << zeroBinWidthPara  << "     " << deltaKPara*(.5) << endl;
		for (int n = 0; n < kParaBins-1; n++){
			outfile4 << deltaKPara*(n+.5) << "     " << deltaKPara*(n+1.5) << endl;
		}
		outfile5 << zeroBinWidthPerp << "     " << deltaKPerp*(.5) << endl;
		for (int n = 0; n < kPerpBins-1; n++){
			outfile5 << deltaKPerp*(n+.5) << "     " << deltaKPerp*(n+1.5) << endl;
		}
		outfile6 << zeroSphericalBinWidth  << "     " << deltaKSpherical*(.5) << endl;
		for (int n = 0; n < kPerpBins-1; n++){
			outfile6 << deltaKSpherical*(n+.5) << "     " << deltaKSpherical*(n+1.5) << endl;
		}
		outfile4.close();
		outfile5.close();
		outfile6.close();
	}
}

void saveRealQhat(CVector qHat, vector<double>& sphericalPower, int field){
	//Cylidrical
	ofstream PSoutfile;
	stringstream ss1;
	ss1 << qHatsDirectory;
	if (useCrossPower) ss1 << "Cross_";
	ss1 << "PSE_field" << field << ".dat";
	PSoutfile.open((ss1.str()).c_str(), ios::trunc);
	for (int j = 0; j < kParaBins*kPerpBins; j++) PSoutfile << qHat.real[j] << endl;
	PSoutfile.close();
	
	//Spherical
	ofstream sphericalPSoutfile;
	stringstream ss2;
	ss2 << sphericalQHatsDirectory;
	if (useCrossPower) ss2 << "Cross_";
	ss2 << "Spherical_PSE_field" << field << ".dat";
	sphericalPSoutfile.open((ss2.str()).c_str(), ios::trunc);
	for (int j = 0; j < kSphereBins; j++) sphericalPSoutfile << sphericalPower[j] << endl;
	sphericalPSoutfile.close();
}

void saveBandPower(CVector& P, int& mostRecentIteration, int& thisInstanceCounter, int thisInstance, vector<double>& sphericalPower){
	/*string iterationCounterFilename = "../QuadraticEstimators/iterationCounter.txt";
	int currentCount;
	while (true){
		fstream iterationInfile(iterationCounterFilename.c_str(),fstream::in);
		iterationInfile >> currentCount;
		iterationInfile.close();
		if (currentCount >= mostRecentIteration){
			mostRecentIteration = currentCount;
			break;
		}
	}	
	ofstream iterationOutfile;
	iterationOutfile.open(iterationCounterFilename.c_str(),ios::trunc);
	iterationOutfile << currentCount + 1;
	iterationOutfile.close();*/

	//Count the number of PSE files so far
	FILE *in;
  	char buff[512];
  	string findString = "find " + qHatsDirectory  + "[0-9]* -nowarn -type  f | wc -l";
  	if(!(in = popen(findString.c_str(), "r"))){
  		return;
  	}
	
	int nFiles;
  	while(fgets(buff, sizeof(buff), in)!=NULL){
  	  	if (isdigit(buff[0])){ //handles the case where there are no such files, and thus assumes that mostRecenterIteration = 0
  			nFiles = atoi(buff);
	  	} else {
	  		nFiles = 0;
		}
		break;
  	}
	pclose(in);
	mostRecentIteration = nFiles;
	cout << "Most Recent Interation = " << mostRecentIteration << endl;


	ofstream PSoutfile;
	stringstream ss;
	ss << qHatsDirectory << thisInstance << "/pse_" << thisInstanceCounter << ".dat";
	cout << ss.str() << endl;
	PSoutfile.open((ss.str()).c_str(), ios::trunc);
	for (int j = 0; j < kParaBins*kPerpBins; j++) PSoutfile << P.real[j] << endl;
	PSoutfile.close();
	
	ofstream sphericalPSoutfile;
	stringstream ss2;
	ss2 << sphericalQHatsDirectory << thisInstance << "/spherical_pse_" << thisInstanceCounter << ".dat";
	cout << ss2.str() << endl;
	sphericalPSoutfile.open((ss2.str()).c_str(), ios::trunc);
	for (int j = 0; j < kSphereBins; j++) sphericalPSoutfile << sphericalPower[j] << endl;
	sphericalPSoutfile.close();
	thisInstanceCounter++;
}


int main(int argc, char *argv[]){
	//------------------Set-up----------------------------//
	s = new Specs();
	srand(time(NULL)+42*atoi(argv[1]));
	loadSpecs("../Specifications/dataSpecs.txt", "../cubeParameters.txt", "../Specifications/powerSpectrumSpecs.txt","../Specifications/powerSpectrumAdvancedOptions.txt");
	int mostRecentIteration = 0;
	int thisInstantceCounter = 0;
	int thisInstance = atoi(argv[1]);
	
	//------------------Load in Covariance Matricies----------------------------//
	cout << endl << "----------------------------------------" << endl << endl;
	cout << "Now loading in covariance matrix information..." << endl << endl;
	CVector mask = loadRealData(dataCubeMaskFilename);
	ResolvedCovariance R;
	loadR(R, "../CovarianceMatrices/R/RSources.dat", "../CovarianceMatrices/R/RToeplitz/RT", ".dat", "../CovarianceMatrices/R/RToeplitz/RTEigen", 
		".dat", "../CovarianceMatrices/R/RToeplitz/RTEigenMasked", ".dat","../CovarianceMatrices/R/RDcross.dat", "../CovarianceMatrices/R/RDmean.dat");
	UnresolvedCovariance U;
	loadU(U, "../CovarianceMatrices/U/Uinfo.dat", "../CovarianceMatrices/U/UX.dat", "../CovarianceMatrices/U/UY.dat",  "../CovarianceMatrices/U/UF.dat", "../CovarianceMatrices/U/UDcross.dat", "../CovarianceMatrices/U/UDmean.dat",
		"../CovarianceMatrices/U/UxEigen.dat", "../CovarianceMatrices/U/UyEigen.dat", "../CovarianceMatrices/U/UfEigen.dat", "../CovarianceMatrices/U/UfdEigen.dat", "../CovarianceMatrices/U/UfEigenMasked.dat", "../CovarianceMatrices/U/UfdEigenMasked.dat");
	GalacticCovariance G;
	loadG(G,"../CovarianceMatrices/G/Ginfo.dat", "../CovarianceMatrices/G/GX.dat", "../CovarianceMatrices/G/GY.dat",  "../CovarianceMatrices/G/GF.dat", "../CovarianceMatrices/G/GDcross.dat", "../CovarianceMatrices/G/GDmean.dat", 
		"../CovarianceMatrices/G/GxEigen.dat", "../CovarianceMatrices/G/GyEigen.dat", "../CovarianceMatrices/G/GfEigen.dat", "../CovarianceMatrices/G/GfdEigen.dat", "../CovarianceMatrices/G/GfEigenMasked.dat", "../CovarianceMatrices/G/GfdEigenMasked.dat"); 
	NoiseCovariance N;
	loadN(N, "../CovarianceMatrices/N/N.dat",mask);	
	
	CVector wedgeMask = calculateWedgeMask();

	if (useCrossPower){  //Cross power ===> half the time
		for (int n = 0; n < nElements; n++){
			N.effCov.entry[n] *= sqrt(2); // I believe this should have been 2 rather than sqrt(2)
			N.avgCov.entry[n] *= sqrt(2);
			if (N.Cov.entry[n] != -1) N.Cov.entry[n] *= sqrt(2);
			/*N.effCov.entry[n] *= 2; 
			N.avgCov.entry[n] *= 2;
			if (N.Cov.entry[n] != -1) N.Cov.entry[n] *= 2;*/
		} // Now I believe I wasn't supposed to be correcting here at all...
	}
	if (CequalsI){
		for (int n = 0; n < nElements; n++){
			N.effCov.entry[n] = 1;
			N.Cov.entry[n] = 1;
			N.avgCov.entry[n] = 1;
		}
		cout << "Covariance matrix taken to be the identity." << endl;
	}

	if (JEquals1) cout << "J factor taken to be 1." << endl;
	if (MatrixRecon) MatrixReconstruction(R, U, G, N, mask);
	if (CovTesting) CovarianceTesting(R, U, G, N, "C");
	
	CVector dummyCVector(s);
	dummyCVector.nElements = 1;
	
	//------------------Consturct Preconditioner----------------------------//
	Preconditioner precon;
 	if (!turnOffPreconditioner && !CequalsI && !CequalsN && !CequalsNU){
		constructPreconditioner(precon, R, U, G, N, mask);
	} 

	//------------------Main Algorithm----------------------------//
	
	if (argc > 1){
		if (thisInstance < fields){
			saveKBins();
			int field = thisInstance;
			cout << "Now analyzing the real data for field " << field << "..." << endl << endl;
			CVector temps0 = loadFileIntoCVector(field, 0, "temps");
			//CVector uvweights0 = loadFileIntoCVector(field, 0, "uvweights");
			CVector temps1 = loadFileIntoCVector(field, 1, "temps");
			//CVector uvweights1 = loadFileIntoCVector(field, 1, "uvweights");
			if (useCrossPower){
				cout << "Using cross power!" << endl;
				if (removeMean) subtractMean(temps0,R,U,G);
				if (removeMean) subtractMean(temps1,R,U,G);
				projectOutMaskedOrUnobservedModes(temps0,mask,N,wedgeMask);
				projectOutMaskedOrUnobservedModes(temps1,mask,N,wedgeMask);
				CVector yReal0 = temps0;
				if (!CequalsI) yReal0 = conjugateGradientAlgorithm(temps0, R, U, G, N, precon, mask, wedgeMask);
				CVector yReal1 = temps1;
				if (!CequalsI) yReal1 = conjugateGradientAlgorithm(temps1, R, U, G, N, precon, mask, wedgeMask);
				projectOutMaskedOrUnobservedModes(yReal0,mask,N,wedgeMask);
				projectOutMaskedOrUnobservedModes(yReal1,mask,N,wedgeMask);
				vector<double> realSphericalPower(kSphereBins,0);
				CVector qHatReal = bandPowerSpectrum(yReal0, yReal1, N, true, realSphericalPower);
				saveRealQhat(qHatReal, realSphericalPower, field);
			} else {
				CVector realData = temps0 + temps1;
				for (int n = 0; n < nElements; n++) realData.real[n] /= 2;
				if (removeMean) subtractMean(realData,R,U,G);
				projectOutMaskedOrUnobservedModes(realData,mask,N,wedgeMask);
				CVector yReal = realData;
				if (!CequalsI) yReal = conjugateGradientAlgorithm(realData, R, U, G, N, precon, mask, wedgeMask);
				projectOutMaskedOrUnobservedModes(yReal,mask,N,wedgeMask);
				vector<double> realSphericalPower(kSphereBins,0);
				CVector qHatReal = bandPowerSpectrum(yReal, dummyCVector, N, true, realSphericalPower);
				saveRealQhat(qHatReal, realSphericalPower, field);
			}
			fftw_cleanup();			
		}
	}
	while (mostRecentIteration < FisherMCNum){
		CVector x = generateUniverse(R,U,G,N);
		projectOutMaskedOrUnobservedModes(x,mask,N,wedgeMask);
		CVector y = x;	
		if (!CequalsI) y = conjugateGradientAlgorithm(x, R, U, G, N, precon, mask, wedgeMask);	
		projectOutMaskedOrUnobservedModes(y,mask,N,wedgeMask);
		vector<double> sphericalPower(kSphereBins,0);
		CVector qHat = bandPowerSpectrum(y, dummyCVector, N, true, sphericalPower);
		if (printAllPSE) saveBandPower(qHat,mostRecentIteration,thisInstantceCounter,thisInstance,sphericalPower);
		fftw_cleanup();
	}

	//TODO: go back to removing the mean and add a function to remove the model mean from the data.  Make sure that this still reporduces Adrian's fisher matrix.
	
	return 0;
}
