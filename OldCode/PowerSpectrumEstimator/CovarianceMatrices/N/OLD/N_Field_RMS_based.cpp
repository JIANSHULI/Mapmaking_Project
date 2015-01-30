#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <map>
#include "../../CommonClasses/Specs.h"
#include <time.h>
#include "../../CommonClasses/CVector.h"
#include <sstream>

using namespace std;

//Constants and global variables
const double pi = 3.1415926535897932384626433832795;
const double c = 299792000; //m/s
const double H0 = 70900; //m/s/Mpc
const double OmegaM = .27;
const double OmegaL = .73;
const double f21cm = 1420.41; //MHz
double deltaRedshift = .00001;
int comovingDistIntSteps = 10000;
Specs *s;
double effectiveArea, systemTemperature, observationTime, fLength, xyLength, fStart;
bool useEmpiricalSigma;
string UVWeightsFilename, maskDataCubeFilename, temps1CubeFilename, temps2CubeFilename;
string cubeDirectory = "../../Cubes/";
string dataCubeMaskFilename = "../../Cubes/dataCubeMask.dat";
string dataCubeString1 = "_field_";
string datacubeString2 = "_slice_";
string datacubeString3 = ".dat";
int xBins, yBins, fBins, nElements, nAntennas, fields, centerField;
double fractionOfModesToThrowOut = .0;
bool ignoreTopAndLeftEdge = false;
double throwOutModesWithObservationTimeRationsGreaterThanThis = 10;

void loadSpecs(string dataSpecsFilename, string cubeParametersFilename, string NSpecsFilename){
	//Load in data files
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
	
	//Load relevant specifications into Specs object
	s->xBins = xBins;
	s->yBins = yBins;
	s->fBins = fBins;
	s->fLength = fLength;
	s->fStart = fStart;
	nElements = xBins*yBins*fBins;

	//Load in noise specs
	infile.open(NSpecsFilename.c_str(),fstream::in);
	infile >> dummy >> effectiveArea;
	infile >> dummy >> systemTemperature;
	infile >> dummy >> nAntennas;
	infile >> dummy >> observationTime;
	infile >> dummy >> useEmpiricalSigma;
	infile.close();	
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

void deconvolve(CVector& dataCube, CVector& UVweights){
	CVector UVweightsFT = UVweights.ijk2uvk();
	CVector FTdata = dataCube.ijk2uvk();	
	CVector FTdataDecon = FTdata;
	for (int n = 0; n < nElements; n++){
		if ((pow(UVweightsFT.real[n],2) + pow(UVweightsFT.imag[n],2)) > 1e-12){
			FTdataDecon.real[n] = (FTdata.real[n]*UVweightsFT.real[n] + FTdata.imag[n]*UVweightsFT.imag[n])/(pow(UVweightsFT.real[n],2) + pow(UVweightsFT.imag[n],2));
			FTdataDecon.imag[n] = (-FTdata.real[n]*UVweightsFT.imag[n] + FTdata.imag[n]*UVweightsFT.real[n])/(pow(UVweightsFT.real[n],2) + pow(UVweightsFT.imag[n],2));
		}
	}
	dataCube = FTdataDecon.uvk2ijk();
	for (int n = 0; n < nElements; n++)	dataCube.imag[n] = 0;
}

//For a given slice, loads in all temps and uvweights and deconvolves them all
vector<CVector> loadAllDeconvolvedFields(int slice){
	vector<CVector> deconvolvedTemps;
	for (int f = 0; f < fields; f++){
		CVector tempsHere = loadFileIntoCVector(f, slice, "temps");
		CVector uvweightsHere = loadFileIntoCVector(f, slice, "uvweights");
		deconvolve(tempsHere,uvweightsHere);
		deconvolvedTemps.push_back(tempsHere);
	}
	return deconvolvedTemps;
}

//Calculate variances in uvk 
vector< vector< vector<double> > > calculateMOSinUV(){
	vector< vector< vector<double> > > MeanOfSquares(xBins, vector< vector<double> >(yBins,vector<double>(fBins,0)));
	vector<CVector> deconvolvedTemps0 = loadAllDeconvolvedFields(0);
	vector<CVector> deconvolvedTemps1 = loadAllDeconvolvedFields(1);
	vector<CVector> allFieldsDiff = deconvolvedTemps0;
	for (int f = 0; f < fields; f++){
		allFieldsDiff[f] = (deconvolvedTemps0[f]-deconvolvedTemps1[f]).ijk2uvk();
	}
	
	for (int f = 0; f < fields; f++){
		for (int k = 0; k < fBins; k++){	
			for (int i = 0; i < xBins; i++){
				for (int j = 0; j < yBins; j++){
					MeanOfSquares[i][j][k] += (pow(allFieldsDiff[f].real[i*yBins*fBins + j*fBins + k],2) + pow(allFieldsDiff[f].imag[i*yBins*fBins + j*fBins + k],2)) / (2.0*fields * xBins * yBins);
					//The 2 is for computing a noise covariance of the mean of slices 0 and 1.
				}
			}
		}
	}
	return MeanOfSquares;
}

//Throws out origin and cells that have high variability in their RMS values
void maskCells(vector< vector< vector<double> > >& MeanOfSquares, vector< vector< vector<double> > >& observationTimes){
	int highRatioModesThrownOut = 0;
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			double smallest = 1e100;
			double largest = 0;
			for (int k = 0; k < fBins; k++){
				if (observationTimes[i][j][k] > 0){
					if (MeanOfSquares[i][j][k] > largest) largest = MeanOfSquares[i][j][k];
					if (MeanOfSquares[i][j][k] < smallest) smallest = MeanOfSquares[i][j][k];
				}
			}
			if (largest/smallest > throwOutModesWithObservationTimeRationsGreaterThanThis){
				highRatioModesThrownOut++;
				for (int k = 0; k < fBins; k++) observationTimes[i][j][k] = 0;
				for (int k = 0; k < fBins; k++) MeanOfSquares[i][j][k] = 0;
			}
			if ((i == xBins/2) && (j == yBins/2)){
				for (int k = 0; k < fBins; k++) observationTimes[i][j][k] = 0;
				for (int k = 0; k < fBins; k++) MeanOfSquares[i][j][k] = 0;
			}
		}
	}
	cout << endl << highRatioModesThrownOut << " modes thrown out due to highly inconsistant observation times." << endl << endl;
}

void enforceSymmetry(vector< vector< vector<double> > >& Cov){
	for (int k = 0; k < fBins; k++){
		for (int i = 0; i < xBins; i++){
			for (int j = 0; j < yBins; j++){
				if ((i == 0) && !(j == 0)){
					double covHere = Cov[i][j][k];
					double covThere = Cov[i][yBins - j][k];
					if (covHere == 0 || covThere == 0){
						Cov[i][j][k] = 0;
						Cov[i][yBins - j][k] = 0;	
					} else {
						Cov[i][j][k] = (covHere + covThere)/2;	
						Cov[i][yBins - j][k] = (covHere + covThere)/2;	
					}
				} else if ((j == 0) && !(i == 0)){
					double covHere = Cov[i][j][k];
					double covThere = Cov[xBins - i][j][k];
					if (covHere == 0 || covThere == 0){
						Cov[i][j][k] = 0;
						Cov[xBins - i][j][k] = 0;	
					} else {
						Cov[i][j][k] = (covHere + covThere)/2;	
						Cov[xBins - i][j][k] = (covHere + covThere)/2;	
					}
				} else if (!(j == 0 && i == 0)){
					double covHere = Cov[i][j][k];
					double covThere = Cov[xBins - i][yBins - j][k];
					if (covHere == 0 || covThere == 0){
						Cov[i][j][k] = 0;
						Cov[xBins - i][yBins - j][k] = 0;	
					} else {
						Cov[i][j][k] = (covHere + covThere)/2;	
						Cov[xBins - i][yBins - j][k] = (covHere + covThere)/2;	
					}
				}
			}
		}
	}
}

void displayEmpiricalRMS(vector< vector< vector<double> > >& MeanOfSquares){
	int countMaskedChannels = 0;
	for (int k = 0; k < fBins; k++){
		if (MeanOfSquares[4][4][k]==0) countMaskedChannels++;
	}

	double overallRMS = 0;
	for (int k = 0; k < fBins; k++){	
		for (int i = 0; i < xBins; i++){
			for (int j = 0; j < yBins; j++){
				overallRMS += MeanOfSquares[i][j][k] / (fBins - countMaskedChannels) / xBins / yBins;
			}
		}
	}
	cout << "Empirical RMS is " << sqrt(overallRMS) << endl;
}


//The average of the synthesized beams for the two time slices is taken to be the synthesized beam from which we determine relative uv coverage
CVector loadCenterFieldAverageUVweights(){
	CVector uvweights0 = loadFileIntoCVector(centerField, 0, "uvweights");
	CVector uvweights1 = loadFileIntoCVector(centerField, 1, "uvweights");
	CVector uvweightsAvg = uvweights0 + uvweights1;
	for (int n = 0; n < nElements; n++) uvweightsAvg.real[n] /= 2;
	return uvweightsAvg;
}

int countMaskedChannels(string filename){
	int maskedChannels = 0;
	fstream infile(filename.c_str(),fstream::in);
	int masked = 0;
	for (int k = 0; k < fBins; k++){
		infile >> masked;
		maskedChannels += masked;
	}
	infile.close();
	return maskedChannels;
}

void printObsTimesToFile(vector< vector< vector<double> > >& observationTimes, int nMaskedChannels){
	vector< vector<double> > avgObsTimes(xBins, vector<double>(yBins,0));
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			for (int k = 0; k < fBins; k++){
				avgObsTimes[i][j] += observationTimes[i][j][k] / (fBins - nMaskedChannels); //averaged observation times
			}
		}
	}
	
	ofstream outfile;
	string ObservationTimesFilename = "obsTimes.dat";
	outfile.open(ObservationTimesFilename.c_str(), ios::trunc);	
	for (int j = 0; j < yBins; j++){
		for (int i = 0; i < xBins; i++) outfile << avgObsTimes[i][j] << " ";
		outfile << endl;
	}
	outfile.close();		
}

vector< vector< vector<double> > > calculateObservationTimes(CVector& uvweights, int nMaskedChannels){
	//Convert to uvk basis (Fourier perpendicular, real parallel)
	CVector FTWeights = uvweights.ijk2uvk();
	for (int n = 0; n < nElements; n++){
		FTWeights.real[n] = sqrt(pow(FTWeights.real[n],2) + pow(FTWeights.imag[n],2));
		FTWeights.imag[n] = 0.0;
	}
	vector< vector< vector<double> > > obsTimes(xBins, vector< vector<double> >(yBins,vector<double>(fBins,0)));
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			for (int k = 0; k < fBins; k++){
				obsTimes[i][j][k] = FTWeights.real[i*yBins*fBins + j*fBins + k];
			}
		}
	}

	vector< vector<double> > avgObsTimes(xBins, vector<double>(yBins,0));
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			for (int k = 0; k < fBins; k++){
				avgObsTimes[i][j] += FTWeights.real[i*yBins*fBins + j*fBins + k] / (fBins - nMaskedChannels); //averaged observation times
				//if (i == 14 && j == 14) cout << FTWeights.real[i*yBins*fBins + j*fBins + k] << endl;
			}
		}
	}
	
	//Normalize and convert to seconds
	double obsTimesSum = 0;
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++) obsTimesSum += avgObsTimes[i][j];
	}
	//cout << "obsTimesSum: " << obsTimesSum << endl;
	double converstionToSeconds = observationTime * 3600.0 * nAntennas * (nAntennas - 1) / obsTimesSum;
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			for (int k = 0; k < fBins; k++){
				obsTimes[i][j][k] *= converstionToSeconds;
			}
			avgObsTimes[i][j] *= converstionToSeconds;
			if (ignoreTopAndLeftEdge){
				if (i == 0 || j == 0){
					for (int k = 0; k < fBins; k++){
						obsTimes[i][j][k] = 0;
					}	
					avgObsTimes[i][j] = 0;
				}
			}
		}
	}

	// Set Origin ObsTimes to 0
	for (int k = 0; k < fBins; k++){
		obsTimes[xBins/2][yBins/2][k] = 0;	
	}
	
	return obsTimes;	
}

vector<double> listFrequencies(){
	double comovingDist = 0;
	double zLeft = f21cm/fStart - 1;
	double zRight = zLeft + deltaRedshift;
	while (true){
		comovingDist += c/H0*((1.0/sqrt(OmegaM*pow(1+zLeft,3) + OmegaL) + 4.0/sqrt(OmegaM*pow(1+(zLeft + zRight)/2,3) + OmegaL) + 1.0/sqrt(OmegaM*pow(1+zRight,3)+OmegaL))*deltaRedshift/6);
		if (comovingDist >= fLength) break;
		zLeft = zRight;
		zRight = zLeft + deltaRedshift;
	}
	double fL = f21cm/(zRight + 1);
	double fEnd = fL + (fStart - fL)/fBins;
	vector<double> freqs(fBins,0);
	double deltaF = (fEnd - fStart)/(fBins - 1);
	for (int i = 0; i < fBins; i++) freqs[i] = fStart + i*deltaF;
	return freqs;
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
	return Omega / xBins / yBins;
}

double calculateTMax(vector< vector< vector<double> > >& observationTimes){
	double tMax = 0;
	for (int k = 0; k < fBins; k++){
		for (int j = 0; j < yBins; j++){
			for (int i = 0; i < xBins; i++){
				if (observationTimes[i][j][k] > tMax) tMax = observationTimes[i][j][k];
			}
		}
	}
	return tMax;
}

void enforceObsTimeSymmetry(vector< vector< vector<double> > >& observationTimes){
	for (int k = 0; k < fBins; k++){
		for (int i = 0; i < xBins; i++){
			for (int j = 0; j < yBins; j++){
				if ((i == 0) && !(j == 0)){
					double timeHere = observationTimes[i][j][k];
					double timeThere = observationTimes[i][yBins - j][k];
					if (timeHere == 0 || timeThere == 0){
						observationTimes[i][j][k] = 0;
						observationTimes[i][yBins - j][k] = 0;	
					} else {
						observationTimes[i][j][k] = (timeHere + timeThere)/2;	
						observationTimes[i][yBins - j][k] = (timeHere + timeThere)/2;	
					}
				} else if ((j == 0) && !(i == 0)){
					double timeHere = observationTimes[i][j][k];
					double timeThere = observationTimes[xBins - i][j][k];
					if (timeHere == 0 || timeThere == 0){
						observationTimes[i][j][k] = 0;
						observationTimes[xBins - i][j][k] = 0;	
					} else {
						observationTimes[i][j][k] = (timeHere + timeThere)/2;	
						observationTimes[xBins - i][j][k] = (timeHere + timeThere)/2;	
					}
				} else if (!(j == 0 && i == 0)){
					double timeHere = observationTimes[i][j][k];
					double timeThere = observationTimes[xBins - i][yBins - j][k];
					if (timeHere == 0 || timeThere == 0){
						observationTimes[i][j][k] = 0;
						observationTimes[xBins - i][yBins - j][k] = 0;	
					} else {
						observationTimes[i][j][k] = (timeHere + timeThere)/2;	
						observationTimes[xBins - i][yBins - j][k] = (timeHere + timeThere)/2;	
					}
				}
			}
		}
	}
}

void updateN(int n, double deltaF, double omegaPix, double f, vector< vector< vector<double> > >& observationTimes, double tMax, vector< vector < vector<double> > >& noiseCovarianceMatrixDiagonal){
	double DeltaX = xyLength / xBins;
	double DeltaY = xyLength / yBins;
	double DeltaKx = 2*pi / xyLength;
	double DeltaKy = 2*pi / xyLength;
	vector<double> kX(xBins,0);
	vector<double> kY(yBins,0);
	for (int i = -xBins/2; i < xBins/2; i++) kX[i + xBins/2] = DeltaKx*(i);
	for (int j = -yBins/2; j < yBins/2; j++) kY[j + yBins/2] = DeltaKy*(j);
	
	double prefactor = (pow(c/1e6/f,4) * pow(systemTemperature,2)) / (pow(effectiveArea,2) * pow(omegaPix,2) * xBins * yBins * (deltaF*1e6));
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			double jx = 1;
			double jy = 1;
			if (kX[i] != 0) jx = sin(kX[i]*DeltaX/2) / (kX[i]*DeltaX/2);
			if (kY[j] != 0) jy = sin(kY[j]*DeltaY/2) / (kY[j]*DeltaY/2);
			if (observationTimes[i][j][n] == 0 || (i == xBins/2 && j == xBins/2)){
				noiseCovarianceMatrixDiagonal[i][j][n] = -1;
			} else {
				noiseCovarianceMatrixDiagonal[i][j][n] =  pow(jx*jy,2) * prefactor / observationTimes[i][j][n];
			}
		}
	}
}


void printTheorySigma(vector< vector < vector<double> > >& noiseCovarianceMatrixDiagonal, vector< vector< vector<double> > >& observationTimes){
	double oldSigmaSquared = 0;
	int maskedElementsCounter = 0;
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			for (int k = 0; k < fBins; k++){
				if (noiseCovarianceMatrixDiagonal[i][j][k] > 0) oldSigmaSquared += noiseCovarianceMatrixDiagonal[i][j][k];
				else maskedElementsCounter++;
			}
		}
	}

	oldSigmaSquared /= (xBins*yBins*fBins - maskedElementsCounter);
	cout << "Theoretical Sigma: " << sqrt(oldSigmaSquared) << endl;
}


void printNtoFile(vector< vector < vector<double> > >& noiseCovarianceMatrixDiagonal){
	string NOutputFilename = "N.dat";
	ofstream outfile;
	outfile.precision(30);
	outfile.open(NOutputFilename.c_str(), ios::trunc);	
	for (int u = 0; u < xBins; u++){
		for (int v = 0; v < yBins; v++){
			for (int k = 0; k < fBins; k++){
				if (noiseCovarianceMatrixDiagonal[u][v][k] > 0) outfile << noiseCovarianceMatrixDiagonal[u][v][k] << endl;
				else outfile << -1 << endl;
			}
		}
	}
	outfile.close();
}

int main(){
	//Load in parameters and data
	cout << endl << "Now calculating the noise covariance matrix..." << endl;
	s = new Specs();
	loadSpecs("../../Specifications/dataSpecs.txt","../../cubeParameters.txt","../../Specifications/NSpecs.txt");
	
	//Calulate Emprical Variances
	vector< vector< vector<double> > > MeanOfSquaresInUVK = calculateMOSinUV();

	//Load in Observation Time information
	int nMaskedChannels = countMaskedChannels(dataCubeMaskFilename);
	CVector centerFieldAverageUVweights = loadCenterFieldAverageUVweights();
	vector< vector< vector<double> > > observationTimes = calculateObservationTimes(centerFieldAverageUVweights,nMaskedChannels);

	//Adjust and Save Empirical Variances
	maskCells(MeanOfSquaresInUVK, observationTimes);
	enforceSymmetry(MeanOfSquaresInUVK);
	displayEmpiricalRMS(MeanOfSquaresInUVK);
	printNtoFile(MeanOfSquaresInUVK);

	//Adjust and save observation times
	enforceSymmetry(observationTimes);
	printObsTimesToFile(observationTimes,nMaskedChannels);

	//Calculate wavelengths and geometric factors
	vector<double> freqs = listFrequencies();
	double deltaF = freqs[0] - freqs[1];
	double bandwidth = fBins * deltaF; 
	double omegaPix = calculateOmegaPix(freqs);
	
	//Calculate N to compare RMS
	vector< vector < vector<double> > > noiseCovarianceMatrixDiagonal(xBins, vector< vector<double> >(yBins, vector<double>(fBins,0)));
	double tMax = calculateTMax(observationTimes);
	for (int n = 0; n < fBins; n++) updateN(n,deltaF,omegaPix,freqs[n],observationTimes,tMax,noiseCovarianceMatrixDiagonal);
	printTheorySigma(noiseCovarianceMatrixDiagonal, observationTimes);
	cout << "Noise covariance calculation complete." << endl << endl;
	return 0;
}

