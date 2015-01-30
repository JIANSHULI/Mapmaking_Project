#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "../../CommonClasses/Specs.h"
#include <string>
#include <sstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <time.h>

using namespace std;

const double pi = 3.1415926535897932384626433832795;
const double c = 299792000; //m/s
const double H0 = 67770; //m/s/Mpc
const double OmegaM = .3086;
const double OmegaL = .6914;
const double f21cm = 1420.41; //MHz

double deltaRedshift = .00001;
int comovingDistIntSteps = 1000;
double highestPossibleRandomFlux = 1000; //Jy TODO: Verify the validity of this number

string ROutputFilenamePrefix = "./RToeplitz/RT";
string ROutputFilenameSuffix =  ".dat";
string REigenFilenamePrefix = "./RToeplitz/RTEigen";
string REigenMaskedFilenamePrefix = "./RToeplitz/RTEigenMasked";
string REigenFilenameSuffix = ".dat";
string RDmeanOutputFilename = "RDmean.dat";
string RDcrossOutputFilename = "RDcross.dat";
string RSourcesFilename = "RSources.dat";

string maskDataCubeFilename = "../../Cubes/dataCubeMask.dat";
double fluxCut, referenceFrequency, spectralIndexMean, spectralIndexStd, spectralIndexPercentError, fluxPercentError, eigenvalueCut;
double xyLength, fStart, fLength;
int xBins, yBins, fBins, nElements;

void loadSpecs(string cubeParametersFilename, string RSpecsFilename){
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
	
	
	//Load in R specs
	infile.open(RSpecsFilename.c_str(),fstream::in);
	infile >> dummy >> fluxCut;
	infile >> dummy >> referenceFrequency;
	infile >> dummy >> spectralIndexMean;
	infile >> dummy >> spectralIndexStd;
	infile >> dummy >> spectralIndexPercentError;
	infile >> dummy >> fluxPercentError;
	infile >> dummy >> eigenvalueCut;
		//cout << fluxCut << " " << referenceFrequency << " " << spectralIndexMean  << " "<< spectralIndexStd  << " "<< spectralIndexPercentError  << " "<< eigenvalueCut  << " "<< endl; 

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


int calculateNumSources(vector<double>& freqs, double& omegaPix){	
	double f = freqs[0] + (freqs[fBins - 1] - freqs[0])/2;
	double comovingDist = 0;
	double z = f21cm/f - 1;
	for (int i = 0; i < comovingDistIntSteps; i++){
		double zLeft = z*i/comovingDistIntSteps;
		double zRight = z*(i+1)/comovingDistIntSteps;
		comovingDist += (1.0/sqrt(OmegaM*pow(1+zLeft,3) + OmegaL) + 4.0/sqrt(OmegaM*pow(1+(zLeft + zRight)/2,3) + OmegaL) + 1.0/sqrt(OmegaM*pow(1+zRight,3)+OmegaL))*z/comovingDistIntSteps/6;
	}
	comovingDist *= c/H0;
	double Omega = xyLength * xyLength / comovingDist / comovingDist;
	omegaPix = Omega / xBins / yBins;
	double nSources = 0;
	double intLimit = .880;
	if (fluxCut >= .880) intLimit = fluxCut;
	nSources += 1.92192/pow(intLimit,1.51);
	if (fluxCut < .880) nSources += -4.69333 + 4.26426/pow(fluxCut,.75);
	nSources *= Omega * 1000; //Jy to mJy
	return int(ceil(nSources));
}

vector<double> sourceFluxes(int nSources){
	vector<double> fluxes;
	double maxFlux = 0;
	while (fluxes.size() < nSources){
 		double flux = highestPossibleRandomFlux*rand()/RAND_MAX;
		double rejectionCriterion;
		if (fluxCut <= .880) rejectionCriterion = 4*pow(fluxCut/.880,-1.75)*rand()/RAND_MAX;
		else rejectionCriterion = 4*pow(fluxCut/.880,-2.51)*rand()/RAND_MAX;
		if (flux >= fluxCut && flux <= .880 && rejectionCriterion < 4*pow(flux/.880,-1.75)){
			fluxes.push_back(flux);
			if (flux > maxFlux) maxFlux = flux;
		}
		if (flux >= fluxCut && flux > .880 && rejectionCriterion < 4*pow(flux/.880,-2.51)){
			fluxes.push_back(flux);
			if (flux > maxFlux) maxFlux = flux;
		} 

	}
	cout << "The largest point source flux is " << maxFlux << " Jy." << endl;
	return fluxes;
	
}

vector < vector<int> > sourcePositions(int nSources){
	vector < vector<int> > positions(nSources, vector<int>(2,0));
	vector < vector<bool> > usedPositions(xBins, vector<bool>(yBins,false));
	int positionsAssigned = 0;
	while (positionsAssigned < nSources){
		int x = rand()%xBins;
		int y = rand()%yBins;
		if (!usedPositions[x][y]){
			positions[positionsAssigned][0] = x;
			positions[positionsAssigned][1] = y;
			usedPositions[x][y] = true;
			positionsAssigned++;
		}
		if (positionsAssigned == xBins*yBins) {
			cout << endl << "Warning! Every pixel has been assigned a resolved point source. Flux cut is too small." << endl << endl;
			break;
		}
	}
	return positions;
}

vector<double> sourceSpectralIndices(int nSources){
	vector<double> spectralIndices(nSources,0);
	double shallowest = 10000000;
	double steepest = -10000000;
	for (int n = 0; n < nSources; n++){
		double r1 = 1.0*rand()/RAND_MAX;
		double r2 = 1.0*rand()/RAND_MAX;
		spectralIndices[n] = spectralIndexStd * sqrt(-2.0*log(r1))*cos(2.0*pi*r2) + spectralIndexMean;
		if (spectralIndices[n] < shallowest) shallowest = spectralIndices[n]; 
		if (spectralIndices[n] > steepest) steepest = spectralIndices[n];
	}
	cout << "The steepest spectral index is " << steepest << " and the shallowest is " << shallowest << "." << endl;
	return spectralIndices;
}


void printSourcesToFile(int nSources, double omegaPix, vector<double>& fluxes, vector < vector<int> >& positions, vector<double>& spectralIndices, vector<double>& freqs){
	ofstream outfile;
	outfile.precision(30);
	outfile.open(RSourcesFilename.c_str(), ios::trunc);	
	outfile << min(nSources,xBins*yBins) << endl << spectralIndexPercentError << endl << fluxPercentError << endl << referenceFrequency << endl << omegaPix << endl;
	for (int n = 0; n < min(nSources,xBins*yBins); n++){
		outfile << positions[n][0] << "    " << positions[n][1] << "    " << fluxes[n] << "    " << spectralIndices[n] << endl;
	}
	for (int n = 0; n < fBins; n++){
		outfile << freqs[n] << endl;
	}
	outfile.close();
	
	ofstream outfile2;
	string positionsFilename = "RSourcePositions.dat";
	outfile2.open(positionsFilename.c_str(), ios::trunc);
	for (int n = 0; n < min(nSources,xBins*yBins); n++){
		outfile2 << positions[n][0] << "    " << positions[n][1] << "    " << endl;
	}
	outfile2.close();
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

void printToeplitzMatrixToFile(vector<double> toeplitz, string outFilename){
	ofstream outfile;
	outfile.precision(30);
	outfile.open(outFilename.c_str(), ios::trunc);	
	for (int i = 0; i < fBins; i++)	outfile << toeplitz[i] << endl;
	outfile.close();
}

void findAndPrintToeplitzEigenSystem(vector<double> toeplitz, vector<double> RDmean, vector<double> RDcross, string outFilename, vector<int>& mask, bool useMask){
	ofstream outfile;
	outfile.precision(30);
	outfile.open(outFilename.c_str(), ios::trunc);
	
	int length = toeplitz.size();
	double matrix[length*length];
	for (int i = 0; i < length; i++){
		for (int j = 0; j < length; j++){
			if (!useMask || mask[i]==0 && mask[j]==0){
				matrix[i*length + j] = RDcross[i]*toeplitz[abs(i - j)]*RDcross[j] - RDmean[i]*RDmean[j];
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
	time_t start, end;
	time(&start);	
	cout << endl << "Now computing resolved point sources covariance matrix componenets..." << endl;
	loadSpecs("../../cubeParameters.txt","../../Specifications/Rspecs.txt");
	srand(time(NULL));
	
	vector<int> mask = loadMask(maskDataCubeFilename);
	vector<double> freqs = listFrequencies();
	double deltaF = freqs[0] - freqs[1];
	cout << "deltaF: " << deltaF << endl;
	double centerFreq = freqs[fBins/2 - 1]; //This is because I'm going reverse the frequency data before FFTing.
	double omegaPix;
	int nSources = calculateNumSources(freqs, omegaPix);
	cout << nSources << " resolved point sources added. " << 100.0 * nSources / xBins / yBins  << "% of pixels have resolved point sources in them." << endl;
	vector<double> fluxes = sourceFluxes(nSources);
	vector < vector<int> > positions = sourcePositions(nSources);
	vector<double> spectralIndices = sourceSpectralIndices(nSources);
	printSourcesToFile(nSources, omegaPix, fluxes, positions, spectralIndices, freqs);
	
	vector< vector < vector<double> > > RPSCMD_Dmean(xBins, vector< vector<double> >(yBins, vector<double>(fBins,0)));
	vector< vector < vector<double> > > RPSCMD_Dcross(xBins, vector< vector<double> >(yBins, vector<double>(fBins,0)));	
	int percent = 0;
	int nRPS = min(nSources,xBins*yBins);
	
	for (int n = 0; n < nRPS; n++){
		if (100.0*n/nRPS > percent) {
			percent = int(ceil(100.0*n/(nRPS-1)));
			cout << " " << percent << "% complete. \r" << std::flush;
		}
		vector<double> RPSCMT_RcrossToeplitzWhitened(fBins,0);
		int x = positions[n][0];
		int y = positions[n][1];
		double spectralIndexError = spectralIndexPercentError * spectralIndices[n];
		double fluxError = fluxPercentError * fluxes[n];
		for (int m = 0; m < fBins; m++){
			double eta = freqs[m]/referenceFrequency;
			RPSCMD_Dmean[x][y][m] = 1.4e-3 * fluxes[n] * 1/omegaPix * pow(eta,-2-spectralIndices[n]) * exp(pow(spectralIndexError,2)/2 * pow(log(eta),2));
			RPSCMD_Dcross[x][y][m] = 1.4e-3 * sqrt(pow(fluxes[n],2) + pow(fluxError,2)) * 1/omegaPix * pow(eta,-2-spectralIndices[n]) * exp(pow(spectralIndexError,2) * pow(log(eta),2));
			RPSCMT_RcrossToeplitzWhitened[m] = exp(-pow(spectralIndexError/freqs[0],2)/2 * pow(freqs[m]-freqs[0],2));
		}
		stringstream ss1, ss2, ss3;
		ss1 << ROutputFilenamePrefix << n << ROutputFilenameSuffix;
		ss2 << REigenFilenamePrefix << n << REigenFilenameSuffix;
		ss3 << REigenMaskedFilenamePrefix << n << REigenFilenameSuffix;
		printToeplitzMatrixToFile(RPSCMT_RcrossToeplitzWhitened, ss1.str());
		findAndPrintToeplitzEigenSystem(RPSCMT_RcrossToeplitzWhitened, RPSCMD_Dmean[x][y], RPSCMD_Dcross[x][y],ss2.str(),mask,false);
		findAndPrintToeplitzEigenSystem(RPSCMT_RcrossToeplitzWhitened, RPSCMD_Dmean[x][y], RPSCMD_Dcross[x][y],ss3.str(),mask,true);			
	}
	
	printRDMatrixToFile(RPSCMD_Dmean, RDmeanOutputFilename);
	printRDMatrixToFile(RPSCMD_Dcross, RDcrossOutputFilename);
	time(&end);	
	cout << endl << "Total R Calcuation Time: " << difftime(end,start) << " seconds." << endl;
	cout << endl << "RESOLVED POINT SOURCES: DONE." << endl << endl;

	
	return 0;
}
