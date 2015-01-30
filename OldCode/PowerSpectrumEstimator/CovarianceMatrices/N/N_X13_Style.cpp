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
double effectiveArea, systemTemperature, observationTime, fLength, xyLength, fStart, baseRMSonInnerPixels, throwOutModesWithObservationTimeRationsGreaterThanThis;
bool useEmpiricalSigma;
string UVWeightsFilename, maskDataCubeFilename, temps1CubeFilename, temps2CubeFilename;
string cubeDirectory = "../../Cubes/";
string dataCubeMaskFilename = "../../Cubes/dataCubeMask.dat";
string dataCubeString1 = "_field_";
string datacubeString2 = "_slice_";
string datacubeString3 = ".dat";
int xBins, yBins, fBins, nElements, nAntennas, fields, centerField;
bool ignoreTopAndLeftEdge = false;

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
	infile >> dummy >> baseRMSonInnerPixels;
	infile >> dummy >> throwOutModesWithObservationTimeRationsGreaterThanThis;

	infile.close();	
}

//Loads a file in my format into a CVector
CVector loadFileIntoCVector(int field, int slice, string type){
	stringstream ss;
	ss << cubeDirectory << type << dataCubeString1 << field << datacubeString2 << slice << datacubeString3;
	fstream infile((ss.str()).c_str(),fstream::in);
	cout << "Loading " << ss.str() << endl;
	CVector loadedVector(s);
	double value;
	for (int n = 0; n < nElements; n++){
		infile >> value;
		loadedVector.real[n] = value;
	}
	return loadedVector;
}

//The average of the synthesized beams for the two time slices is taken to be the synthesized beam from which we determine relative uv coverage
CVector loadCenterFieldAverageUVweights(){
	CVector uvweights0 = loadFileIntoCVector(centerField, 0, "psf");
	CVector uvweights1 = loadFileIntoCVector(centerField, 1, "psf");
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
	cout << "obsTimesSum: " << obsTimesSum << endl;
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

void deconvolve(CVector& dataCube, CVector& UVweights){
	cout << "Deconvolution turned off." << endl;
	return;
	/*CVector UVweightsFT = UVweights.ijk2uvk();
	CVector FTdata = dataCube.ijk2uvk();	
	CVector FTdataDecon = FTdata;
	for (int n = 0; n < nElements; n++){
		if ((pow(UVweightsFT.real[n],2) + pow(UVweightsFT.imag[n],2)) > 1e-12){
			FTdataDecon.real[n] = (FTdata.real[n]*UVweightsFT.real[n] + FTdata.imag[n]*UVweightsFT.imag[n])/(pow(UVweightsFT.real[n],2) + pow(UVweightsFT.imag[n],2));
			FTdataDecon.imag[n] = (-FTdata.real[n]*UVweightsFT.imag[n] + FTdata.imag[n]*UVweightsFT.real[n])/(pow(UVweightsFT.real[n],2) + pow(UVweightsFT.imag[n],2));
		}
	}
	dataCube = FTdataDecon.uvk2ijk();
	for (int n = 0; n < nElements; n++)	dataCube.imag[n] = 0;*/
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

void maskVeryNoisyCells(vector< vector< vector<double> > >& observationTimes){
	//Now throw out modes where the ratio of the most observeration to the least is very large:
	int highRatioModesThrownOut = 0;
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			double smallestObs = 1e100;
			double largestObs = 0;
			for (int k = 0; k < fBins; k++){
				if (observationTimes[i][j][k] > 0){
					if (observationTimes[i][j][k] > largestObs) largestObs = observationTimes[i][j][k];
					if (observationTimes[i][j][k] < smallestObs) smallestObs = observationTimes[i][j][k];
				}
			}
			double dist = sqrt(1.0*pow((i-xBins/2),2) + 1.0*pow((j-yBins/2),2));
			if (largestObs/smallestObs > throwOutModesWithObservationTimeRationsGreaterThanThis && dist > baseRMSonInnerPixels){
				highRatioModesThrownOut++;
				for (int k = 0; k < fBins; k++) observationTimes[i][j][k] = 0;
			}
		}
	}
	cout << highRatioModesThrownOut << " modes thrown out due to highly inconsistant observation times (ratio > " << throwOutModesWithObservationTimeRationsGreaterThanThis << ")." << endl << endl;
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

double calculateEmpiricalSigma(vector< vector< vector<double> > >& observationTimes){
	vector<CVector> deconvolvedTemps0 = loadAllDeconvolvedFields(0);
	vector<CVector> deconvolvedTemps1 = loadAllDeconvolvedFields(1);
	
	for (int f = 0; f < fields; f++){
		deconvolvedTemps0[f] = deconvolvedTemps0[f].ijk2uvk();
		deconvolvedTemps1[f] = deconvolvedTemps1[f].ijk2uvk();
	}

	int countOfUnmaskedBins = 0;
	for (int f = 0; f < fields; f++){
		for (int i = 0; i < xBins; i++){
			for (int j = 0; j < yBins; j++){
				for (int k = 0; k < fBins; k++){
					double dist = sqrt(1.0*pow((i-xBins/2),2) + 1.0*pow((j-yBins/2),2));
					if (observationTimes[i][j][k] == 0 || dist > baseRMSonInnerPixels){
						deconvolvedTemps0[f].real[i*yBins*fBins + j*fBins + k] = 0;
						deconvolvedTemps0[f].imag[i*yBins*fBins + j*fBins + k] = 0;
						deconvolvedTemps1[f].real[i*yBins*fBins + j*fBins + k] = 0;
						deconvolvedTemps1[f].imag[i*yBins*fBins + j*fBins + k] = 0;
					} else {
						countOfUnmaskedBins++;
					}
				}
			}
		}
	}
	for (int f = 0; f < fields; f++){
		deconvolvedTemps0[f] = deconvolvedTemps0[f].uvk2ijk();
		deconvolvedTemps1[f] = deconvolvedTemps1[f].uvk2ijk();
	}

	double totalEmpiricalSigmaSquared = 0;
	//cout << endl << "We have filtered out the top " << 100*fractionOfModesToThrowOut << "% noisiest modes and declared them 'unobserved'" << endl;
	for (int f = 0; f < fields; f++){
		double empricalSigmaSquared = 0;
		for (int i = 0; i < xBins; i++){
			for (int j = 0; j < yBins; j++){
				for (int k = 0; k < fBins; k++){
					double temp0 = deconvolvedTemps0[f].real[i*yBins*fBins + j*fBins + k];
					double temp1 = deconvolvedTemps1[f].real[i*yBins*fBins + j*fBins + k];
					empricalSigmaSquared += pow(temp0-temp1,2);
				}
			}
		}
		empricalSigmaSquared /= (countOfUnmaskedBins/fields);

		empricalSigmaSquared /= 2;
		cout << "Emprical RMS for Field " << f << " is " << sqrt(empricalSigmaSquared) << endl;
		totalEmpiricalSigmaSquared += empricalSigmaSquared / fields;
	}
	cout << "Overall RMS: " << sqrt(totalEmpiricalSigmaSquared) << endl << endl;

	return sqrt(totalEmpiricalSigmaSquared);
}

void enforceSigma(vector< vector < vector<double> > >& noiseCovarianceMatrixDiagonal, vector< vector< vector<double> > >& observationTimes){
	double oldSigmaSquared = 0;
	int maskedElementsCounter = 0;
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			for (int k = 0; k < fBins; k++){
				double dist = sqrt(1.0*pow((i-xBins/2),2) + 1.0*pow((j-yBins/2),2));
				if (noiseCovarianceMatrixDiagonal[i][j][k] > 0 && dist <= baseRMSonInnerPixels) oldSigmaSquared += noiseCovarianceMatrixDiagonal[i][j][k];
				else maskedElementsCounter++;
			}
		}
	}

	oldSigmaSquared /= (xBins*yBins*fBins - maskedElementsCounter);
	cout << "Theoretical Sigma: " << sqrt(oldSigmaSquared) << endl;
	if (useEmpiricalSigma){
		double empiricalSigma = calculateEmpiricalSigma(observationTimes);
		cout << "Empirical Sigma: " << empiricalSigma << endl;
		for (int i = 0; i < xBins; i++){
			for (int j = 0; j < yBins; j++){
				for (int k = 0; k < fBins; k++){
					if (noiseCovarianceMatrixDiagonal[i][j][k] > 0) noiseCovarianceMatrixDiagonal[i][j][k] *= empiricalSigma*empiricalSigma/oldSigmaSquared;
				}
			}
		}
	}
	double largestEV = 0;
	double smallestEV = 1e50;
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			for (int k = 0; k < fBins; k++){
				double here = noiseCovarianceMatrixDiagonal[i][j][k];
				if (here > 0 && here < smallestEV) smallestEV = here;
				if (here > largestEV) largestEV = here;
			}
		}
	}
	cout << "Largest Noise Eigenvalue: " << largestEV << endl;
	cout << "Smallest Noise Eigenvalue: " << smallestEV << endl;
	cout << "Noise Covariance Condition Number: " << largestEV/smallestEV << endl;
}


void printNtoFile(vector< vector < vector<double> > >& noiseCovarianceMatrixDiagonal){
	string NOutputFilename = "N.dat";
	ofstream outfile;
	outfile.precision(30);
	outfile.open(NOutputFilename.c_str(), ios::trunc);	
	for (int u = 0; u < xBins; u++){
		for (int v = 0; v < yBins; v++){
			for (int k = 0; k < fBins; k++){
				outfile << noiseCovarianceMatrixDiagonal[u][v][k] << endl;
			}
		}
	}
	outfile.close();
}

int main(){
	cout << "THE COVARIANCE HERE IS PROBABLY NOT RIGHT, BUT IT'S WHAT I WAS DOING FOR X13" << endl;

	//Load in parameters and data
	cout << endl << "Now calculating the noise covariance matrix..." << endl;
	s = new Specs();
	loadSpecs("../../Specifications/dataSpecs.txt","../../cubeParameters.txt","../../Specifications/NSpecs.txt");
	CVector centerFieldAverageUVweights = loadCenterFieldAverageUVweights();
	
	//Calculate observation times
	int nMaskedChannels = countMaskedChannels(dataCubeMaskFilename);
	vector< vector< vector<double> > > observationTimes = calculateObservationTimes(centerFieldAverageUVweights,nMaskedChannels);
	maskVeryNoisyCells(observationTimes);
	enforceObsTimeSymmetry(observationTimes);
	printObsTimesToFile(observationTimes,nMaskedChannels);

	//Calculate wavelengths and geometric factors
	vector<double> freqs = listFrequencies();
	double deltaF = freqs[0] - freqs[1];
	double bandwidth = fBins * deltaF; 
	double omegaPix = calculateOmegaPix(freqs);
	
	//Calculate N and save
	vector< vector < vector<double> > > noiseCovarianceMatrixDiagonal(xBins, vector< vector<double> >(yBins, vector<double>(fBins,0)));
	double tMax = calculateTMax(observationTimes);
	for (int n = 0; n < fBins; n++) updateN(n,deltaF,omegaPix,freqs[n],observationTimes,tMax,noiseCovarianceMatrixDiagonal);
	enforceSigma(noiseCovarianceMatrixDiagonal, observationTimes);
	printNtoFile(noiseCovarianceMatrixDiagonal);
	cout << "Noise covariance calculation complete." << endl << endl;
	return 0;
}

