#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <map>
#include "../../CommonClasses/Specs.h"
#include <time.h>
#include "../../CommonClasses/CVector.h"
#include <sstream>
#include <gsl/gsl_statistics_double.h>

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
int xBins, yBins, fBins, nElements, nAntennas, fields, centerField, slices;
bool ignoreTopAndLeftEdge = false;

void loadSpecs(string dataSpecsFilename, string cubeParametersFilename, string NSpecsFilename){
	//Load in data files
	fstream infile(dataSpecsFilename.c_str(),fstream::in);
	string dummy;
	for (int n = 0; n < 3; n++) infile >> dummy >> dummy;
	infile >> dummy >> slices;
	infile >> dummy >> fields;
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
	CVector uvweightsAvg = loadFileIntoCVector(centerField, 0, "uvweights");
	for (int l = 1; l < slices; l++) uvweightsAvg = uvweightsAvg + loadFileIntoCVector(centerField, l, "uvweights");
	for (int n = 0; n < nElements; n++) uvweightsAvg.real[n] /= slices;
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

void setAeffAndTsys(vector<double>& freqs){
	cout << endl << "Now adjusting system temperature and effective area to match Aaron" << endl;
	double centerFreq = freqs[fBins/2];
	systemTemperature = 50*pow(centerFreq/150.0,-2.4) + 237*pow(centerFreq/150.0,-2.5);
	cout << "System Temperature set to " << systemTemperature << " K." << endl;
	if (centerFreq < 70){
		cout << "*********************************************************" << endl;
		cout << "WARNING: EFFECTIVE AREA FORMULA NOT VALID FOR F < 70 MHZ." << endl;
		cout << "*********************************************************" << endl;
		centerFreq = 70;
	}
	effectiveArea = 10351.70218945976 - 674.9783701774664*centerFreq + 19.29667765875074*pow(centerFreq,2) - 0.3179099125155949*pow(centerFreq,3) + 0.0033448155237735372*pow(centerFreq,4) - 0.000023505908323372595*pow(centerFreq,5) + 1.118635349336024e-7 *pow(centerFreq,6) - 3.563892774619191e-10 *pow(centerFreq,7) + 7.283615837192128e-13 *pow(centerFreq,8) - 8.633626877485863e-16 *pow(centerFreq,9) + 4.513113181959486e-19 *pow(centerFreq,10);
	cout << "Effective Area set to " << effectiveArea << " m^2." << endl << endl;
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
	prefactor *= 1.0/fields; //This fixes the factor related to the missing 
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			double jx = 1;
			double jy = 1;
			if (kX[i] != 0) jx = sin(kX[i]*DeltaX/2) / (kX[i]*DeltaX/2);
			if (kY[j] != 0) jy = sin(kY[j]*DeltaY/2) / (kY[j]*DeltaY/2);
			if (observationTimes[i][j][n] == 0){
				noiseCovarianceMatrixDiagonal[i][j][n] = -1;
			} else {
				noiseCovarianceMatrixDiagonal[i][j][n] =  pow(jx*jy,2) * prefactor / observationTimes[i][j][n];
			}
		} 
	}
}


vector< vector< vector<double> > > calculateSigmaSqPerUV(){
	vector<CVector> tempsSlices;
	for (int l = 0; l < slices; l++) tempsSlices.push_back(loadFileIntoCVector(centerField, l, "temps"));
	
	vector<CVector> tempsSlicesUV;
	for (int l = 0; l < slices; l++) tempsSlicesUV.push_back(tempsSlices[l].ijk2uvk());

	vector< vector < vector<double> > > UVsigmaSq(xBins, vector< vector<double> >(yBins, vector<double>(fBins,0)));
	for (int k = 0; k < fBins; k++){
		for (int j = 0; j < yBins; j++){
			for (int i = 0; i < xBins; i++){
				double* tempsHereReal = new double[slices];
				double* tempsHereImag = new double[slices];
				for (int l = 0; l < slices; l++){
					tempsHereReal[l] = tempsSlicesUV[l].real[k + j*fBins + i*fBins*yBins];
					tempsHereImag[l] = tempsSlicesUV[l].imag[k + j*fBins + i*fBins*yBins];
				} 
				//for (int l = 0; l < slices; l++) cout << tempsSlicesUV[l].imag[k + j*fBins + i*fBins*yBins] << endl;
				UVsigmaSq[i][j][k] = gsl_stats_variance(tempsHereReal,1,slices)/slices/xBins/yBins + gsl_stats_variance(tempsHereImag,1,slices)/slices/xBins/yBins;
				delete[] tempsHereReal;
				delete[] tempsHereImag;
			}
		}
	}
	return UVsigmaSq;
}



void printNtoFile(vector< vector < vector<double> > >& noiseCovarianceMatrixDiagonal, bool RMS){
	string NOutputFilename = "N.dat";
	if (RMS) NOutputFilename = "N_slices_RMS.dat";

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
	//Load in parameters and data
	cout << endl << "Now calculating the noise covariance matrix..." << endl;
	s = new Specs();
	loadSpecs("../../Specifications/dataSpecs.txt","../../cubeParameters.txt","../../Specifications/NSpecs.txt");
	centerField = 4;
	cout << endl << "centerField assumed to be 4." << endl << endl;
	CVector centerFieldAverageUVweights = loadCenterFieldAverageUVweights();
	
	//Calculate observation times
	int nMaskedChannels = countMaskedChannels(dataCubeMaskFilename);
	vector< vector< vector<double> > > observationTimes = calculateObservationTimes(centerFieldAverageUVweights,nMaskedChannels);
	maskVeryNoisyCells(observationTimes);
	enforceObsTimeSymmetry(observationTimes);
	printObsTimesToFile(observationTimes,nMaskedChannels);

	//Calculate wavelengths and geometric factors
	vector<double> freqs = listFrequencies();
	setAeffAndTsys(freqs);
	double deltaF = freqs[0] - freqs[1];
	double bandwidth = fBins * deltaF; 
	double omegaPix = calculateOmegaPix(freqs);
	
	//Calculate N and save
	vector< vector < vector<double> > > noiseCovarianceMatrixDiagonal(xBins, vector< vector<double> >(yBins, vector<double>(fBins,0)));
	double tMax = calculateTMax(observationTimes);
	for (int n = 0; n < fBins; n++) updateN(n,deltaF,omegaPix,freqs[n],observationTimes,tMax,noiseCovarianceMatrixDiagonal);
	
	//enforceSigma(noiseCovarianceMatrixDiagonal, observationTimes);
	vector< vector< vector<double> > > uvSlicesSigmaSq = calculateSigmaSqPerUV();
	printNtoFile(noiseCovarianceMatrixDiagonal,false);
	printNtoFile(uvSlicesSigmaSq,true);
	cout << "Noise covariance calculation complete." << endl << endl;
	return 0;
}

