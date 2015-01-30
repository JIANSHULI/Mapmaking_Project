#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include "../../CommonClasses/Specs.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h> 
#include "fftw3.h"
#include <time.h>

using namespace std;

//Constants and global variables
const double pi = 3.1415926535897932384626433832795;
const double c = 299792000; //m/s
const double H0 = 67770; //m/s/Mpc
const double OmegaM = .3086;
const double OmegaL = .6914;
const double f21cm = 1420.41; //MHz
double deltaRedshift = .00001;
int comovingDistIntSteps = 10000;
double referenceFrequency = 150;
string USpecsFilename = "Uspecs.txt";
string USourceInfoFilename = "Uinfo.dat";
string UToeplitzXFilename = "UX.dat";
string UToeplitzYFilename = "UY.dat";
string UToeplitzFFilename = "UF.dat";
string UDmeanOutputFilename = "UDmean.dat";
string UDcrossOutputFilename = "UDcross.dat";
string UxEigenFilename = "UxEigen.dat";
string UyEigenFilename = "UyEigen.dat";
string UfEigenFilename = "UfEigen.dat";
string UfdEigenFilename = "UfdEigen.dat";
string UfEigenMaskedFilename = "UfEigenMasked.dat";
string UfdEigenMaskedFilename = "UfdEigenMasked.dat";
Specs *s;
double fLength, xyLength, fStart, fluxCut, spectralIndexMean, spectralIndexStd, clusteringStdInRad, eigenvalueCut;
string maskDataCubeFilename = "../../Cubes/dataCubeMask.dat";
int xBins, yBins, fBins, nElements, nAntennas;

void loadSpecs(string cubeParametersFilename, string USpecsFilename){
	//Load in data cube parameters
	fstream infile(cubeParametersFilename.c_str(),fstream::in);
	string dummy;
	infile >> dummy >> xBins;
	infile >> dummy >> yBins;
	infile >> dummy >> fBins;
	infile >> dummy >> xyLength;
	infile >> dummy >> fLength;
	infile >> dummy >> fStart; 
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
	
	//Load relevant specifications into Specs object
	s->xBins = xBins;
	s->yBins = yBins;
	s->fBins = fBins;
	s->fLength = fLength;
	s->fStart = fStart;
	nElements = xBins*yBins*fBins;
	
	//Load in U specs
	infile.open(USpecsFilename.c_str(),fstream::in);
	infile >> dummy >> fluxCut;
	infile >> dummy >> spectralIndexMean;
	infile >> dummy >> spectralIndexStd;
	infile >> dummy >> clusteringStdInRad;
	infile >> dummy >> eigenvalueCut;
	infile.close();	
}

vector<int> loadMask(string filename){
	vector<int> mask(fBins,0);
	fstream infile(filename.c_str(),fstream::in);
	int value = 0;
	for (int n = 0; n < fBins; n++){
		infile >> value;
		mask[n] = value;
	}
	return mask;
}

vector<double> listFrequencies(){
	double comovingDist = 0;
	double zRight = f21cm/fStart - 1;
	double zLeft = zRight + deltaRedshift;
	while (true){
		comovingDist -= c/H0*((1.0/sqrt(OmegaM*pow(1+zRight,3) + OmegaL) + 4.0/sqrt(OmegaM*pow(1+(zLeft + zRight)/2,3) + OmegaL) + 1.0/sqrt(OmegaM*pow(1+zLeft,3)+OmegaL))*deltaRedshift/6);
		if (comovingDist <= -fLength) break;
		zRight = zLeft;
		zLeft = zRight - deltaRedshift;
	}
	double fL = f21cm/(zLeft + 1);
	double fEnd = fL + (fStart - fL)/fBins;
	vector<double> freqs(fBins,0);
	double deltaF = (fEnd - fStart)/(fBins - 1);
	for (int i = 0; i < fBins; i++){
		freqs[i] = fStart + (fBins - 1 - i)*deltaF;
	} 
	return freqs;
}


vector<double> listXPositions(){
	vector<double> positions(xBins,0);
	double deltaX = xyLength/xBins;
	for (int i = 0; i < xBins; i++) positions[i] = -deltaX*xBins/2 + i*deltaX + deltaX/2;
	return positions;	
}

vector<double> listYPositions(){
	vector<double> positions(yBins,0);
	double deltaY = xyLength/yBins;
	for (int j = 0; j < yBins; j++) positions[j] = -deltaY*yBins/2 + j*deltaY + deltaY/2;
	return positions;
}

double calculateOmegaPix(vector<double>& freqs){	
	double f = fStart + (freqs[fBins -1] - fStart)/2;		
	double comovingDist = 0;
	double z = f21cm/f - 1;
	for (int i = 0; i < comovingDistIntSteps; i++){
		double zLeft = z*i/comovingDistIntSteps;
		double zRight = z*(i+1)/comovingDistIntSteps;
		comovingDist += (1.0/sqrt(OmegaM*pow(1+zLeft,3) + OmegaL) + 4.0/sqrt(OmegaM*pow(1+(zLeft + zRight)/2,3) + OmegaL) + 1.0/sqrt(OmegaM*pow(1+zRight,3)+OmegaL))*z/comovingDistIntSteps/6;
	}
	comovingDist *= c/H0;
	double Omega = xyLength * xyLength / comovingDist / comovingDist;
	cout << "Pixel angular size: " << Omega / xBins / yBins << " steradians." << endl;
	return Omega / xBins / yBins;
}

double calculateBoxAngularSize(vector<double>& freqs){	
	double f = fStart + (freqs[fBins -1] - fStart)/2;		
	double comovingDist = 0;
	double z = f21cm/f - 1;
	for (int i = 0; i < comovingDistIntSteps; i++){
		double zLeft = z*i/comovingDistIntSteps;
		double zRight = z*(i+1)/comovingDistIntSteps;
		comovingDist += (1.0/sqrt(OmegaM*pow(1+zLeft,3) + OmegaL) + 4.0/sqrt(OmegaM*pow(1+(zLeft + zRight)/2,3) + OmegaL) + 1.0/sqrt(OmegaM*pow(1+zRight,3)+OmegaL))*z/comovingDistIntSteps/6;
	}
	comovingDist *= c/H0;
	return xyLength / comovingDist;
}


double sourceCountIntegral1(){
	//return 8263.68;
	if (fluxCut < .880) return 12792.772 * pow(fluxCut,.25);
	else return 12390.4 + 6073.72549 - 5690.380695 * pow(fluxCut,-.51);
}

double sourceCountIntegral2(){
	//return 17270.2;
	if (fluxCut < .880) return 2558.554395 * pow(fluxCut,1.25);
	else return 2180.710397  - 5563.03673 + 5922.6411316 * pow(fluxCut,.49);
}

void printSourceInfoToFile(double omegaPix, double I1, double I2, vector<double>& freqs){
	ofstream outfile;
	outfile.precision(30);		
	outfile.open(USourceInfoFilename.c_str(), ios::trunc);	
	outfile << spectralIndexMean << endl << spectralIndexStd << endl << referenceFrequency << endl << omegaPix << endl << I1 << endl << I2 << endl;
	for (int n = 0; n < fBins; n++){
		outfile << freqs[n] << endl;
	}
	outfile.close();
}

void printRDMatrixToFile(vector< vector < vector<double> > >& matrixDiagonal, string outFilename){
	ofstream outfile;
	outfile.precision(30);
	outfile.open(outFilename.c_str(), ios::trunc);

	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			for (int k = 0; k < fBins; k++){
				outfile << matrixDiagonal[i][j][k] << endl;
			}
		}
	}
	outfile.close();
}

void printToeplitzToFile(vector<double>& toeplitz, string outFilename){
	ofstream outfile;
	outfile.precision(30);
	outfile.open(outFilename.c_str(), ios::trunc);	
	for (int i = 0; i < toeplitz.size(); i++) outfile << toeplitz[i] << endl;
	outfile.close();
}

bool positiveDefCircEmbedding(vector<double>& toeplitz){
	int N = toeplitz.size()*2-1;
	fftw_complex *circIn, *circOut;
	circIn = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	circOut = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	for (int n = 0; n < toeplitz.size(); n++) circIn[n][0] = toeplitz[n];
	for (int n = 0; n < toeplitz.size() - 1; n++) circIn[n + toeplitz.size()][0] = toeplitz[toeplitz.size() - 1 - n];
	for (int n = 0; n < N; n++) circIn[n][1] = 0;
	fftw_plan circFFT = fftw_plan_dft_1d(N, circIn, circOut, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(circFFT);
	for (int n = 0; n < N; n++){
		if (circOut[n][0] < -1e-10) return false;
	}
	return true;	
}
	

void findToeplitzEigendecomp(vector<double>& toeplitz, vector<double>& diag, string outFilename, vector<int>& mask, bool useMask){
	ofstream outfile;
	outfile.precision(30);
	outfile.open(outFilename.c_str(), ios::trunc);
	
	int length = toeplitz.size();
	double matrix[length*length];
	for (int i = 0; i < length; i++){
		for (int j = 0; j < length; j++){
			if (!useMask || mask[i]==0 && mask[j]==0){
				matrix[i*length + j] = toeplitz[abs(i - j)]*diag[i]*diag[j];
			} else {
				matrix[i*length + j] = 0;
			}
		}
	}
	
	gsl_matrix_view m = gsl_matrix_view_array (matrix, length, length);
	gsl_vector *eval = gsl_vector_alloc (length);
	gsl_matrix *evec = gsl_matrix_alloc (length, length);
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (length);
	gsl_eigen_symmv (&m.matrix, eval, evec, w);
	gsl_eigen_symmv_free (w);
	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);
		
	double largestEVal = gsl_vector_get(eval, 0);
	int numEVals = 0;	
	while (gsl_vector_get(eval, numEVals) > largestEVal * eigenvalueCut) {
		numEVals++;
		if (numEVals == length) break;
	}
	outfile << numEVals << endl;
	for (int i = 0; i < numEVals; i++){
		outfile << gsl_vector_get(eval, i) << endl;
		for (int j = 0; j < length; j++){
			outfile << gsl_matrix_get(evec,j,i) << endl;
		}
	}
	
	gsl_vector_free (eval);
	gsl_matrix_free (evec);
	
	outfile.close();
}


int main(){
	//Load in parameters and data
	cout << endl << "Now calculating the unresolved point sources covariance matrix..." << endl;
	s = new Specs();
	loadSpecs("../../cubeParameters.txt","../../Specifications/USpecs.txt");
	vector<double> freqs = listFrequencies();
	vector<double> xPositions = listXPositions();
	vector<double> yPositions = listYPositions();
	double omegaPix = calculateOmegaPix(freqs);
	double boxAngularSize = calculateBoxAngularSize(freqs);
	double clusteringStd = clusteringStdInRad / boxAngularSize * xyLength;
	double I1 = sourceCountIntegral1();
	double I2 = sourceCountIntegral2();
	
	vector< vector < vector<double> > > UPSCMD_Dmean(xBins, vector< vector<double> >(yBins, vector<double>(fBins,0)));
	vector< vector < vector<double> > > UPSCMD_Dcross(xBins, vector< vector<double> >(yBins, vector<double>(fBins,0)));
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){		
			for (int k = 0; k < fBins; k++){
				double eta = freqs[k]/referenceFrequency;
				UPSCMD_Dmean[i][j][k] = 1.4e-3 * pow(eta,-2-spectralIndexMean) * I1 * exp(pow(spectralIndexStd,2)/2 * pow(log(eta),2));
				UPSCMD_Dcross[i][j][k] = 1.4e-3 * pow(I2/omegaPix,.5) * pow(eta,-2-spectralIndexMean) * exp(pow(spectralIndexStd,2) * pow(log(eta),2));
			}
		}
	}
	

	vector<double> UPSCMT_UcrossToeplitzX(xBins,0);
	for (int i = 0; i < xBins; i++)UPSCMT_UcrossToeplitzX[i] = exp(-pow(xPositions[i]-xPositions[0],2)/2/pow(clusteringStd,2));	
	vector<double> UPSCMT_UcrossToeplitzY(yBins,0);
	for (int j = 0; j < yBins; j++) UPSCMT_UcrossToeplitzY[j] = exp(-pow(yPositions[j]-yPositions[0],2)/2/pow(clusteringStd,2));	
	double deltaF = freqs[0] - freqs[1];
	vector<double> UPSCMT_UcrossToeplitzF(fBins,0);
	for (int k = 0; k < fBins; k++) UPSCMT_UcrossToeplitzF[k] = exp(-pow(spectralIndexStd/freqs[0],2)/2 * pow(freqs[k]-freqs[0],2));
	
	vector<int> mask = loadMask(maskDataCubeFilename);
	vector<double> xIdentity(xBins,1);
	vector<double> yIdentity(yBins,1);
	vector<double> fIdentity(fBins,1);
	findToeplitzEigendecomp(UPSCMT_UcrossToeplitzX,xIdentity,UxEigenFilename,mask,false);
	findToeplitzEigendecomp(UPSCMT_UcrossToeplitzY,yIdentity,UyEigenFilename,mask,false);
	findToeplitzEigendecomp(UPSCMT_UcrossToeplitzF,fIdentity,UfEigenFilename,mask,false);
	findToeplitzEigendecomp(UPSCMT_UcrossToeplitzF,UPSCMD_Dcross[0][0],UfdEigenFilename,mask,false);
	findToeplitzEigendecomp(UPSCMT_UcrossToeplitzF,fIdentity,UfEigenMaskedFilename,mask,true);
	findToeplitzEigendecomp(UPSCMT_UcrossToeplitzF,UPSCMD_Dcross[0][0],UfdEigenMaskedFilename,mask,true);
	
	printRDMatrixToFile(UPSCMD_Dmean, UDmeanOutputFilename);
	printRDMatrixToFile(UPSCMD_Dcross, UDcrossOutputFilename);
	printToeplitzToFile(UPSCMT_UcrossToeplitzX, UToeplitzXFilename);
	printToeplitzToFile(UPSCMT_UcrossToeplitzY, UToeplitzYFilename);
	printToeplitzToFile(UPSCMT_UcrossToeplitzF, UToeplitzFFilename);
	printSourceInfoToFile(omegaPix, I1, I2, freqs);
	cout << endl << "UNRESOLVED POINT SOURCES: DONE." << endl << endl;
	

	return 0;
}

