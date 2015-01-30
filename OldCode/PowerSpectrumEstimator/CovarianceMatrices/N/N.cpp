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
const double H0 = 67770; //m/s/Mpc
const double OmegaM = .3086;
const double OmegaL = .6914;
const double f21cm = 1420.41; //MHz
double deltaRedshift = .00001;
int comovingDistIntSteps = 10000;
Specs *s, *sFullField;
double systemTemperatureConstant, systemTemperatureCoeff, systemTemperautreExp, observationTime, throwOutModesThisMuchWorseThanTheBest; 
double fLength, xyLength, fStart, baseRMSonInnerPixels, throwOutModesWithObservationTimeRationsGreaterThanThis, kernelConvolutionLimit;
bool useEmpiricalSigma, approximatelyDownsampleNoiseCovariance;
string UVWeightsFilename, maskDataCubeFilename, temps1CubeFilename, temps2CubeFilename;
string cubeDirectory = "../../Cubes/";
string bigCubeDirectory = "/nfs/eor-09/r1/abrahamn/josh_x16/Cubes/";
string dataCubeMaskFilename = "../../Cubes/dataCubeMask.dat";
string fullFieldDataCubeMaskFilename = "../../Cubes/fullFieldDataCubeMask.dat";
string dataCubeString1 = "_field_";
string datacubeString2 = "_slice_";
string datacubeString3 = ".dat";
int xBins, yBins, fBins, nElements, nAntennas, fields, centerField, fullFieldxBins, fullFieldyBins, downSamplingFactor;
bool ignoreTopAndLeftEdge = false;


void loadSpecs(string dataSpecsFilename, string cubeParametersFilename, string NSpecsFilename){
	//Load in data files
	fstream infile(dataSpecsFilename.c_str(),fstream::in);
	string dummy;
	for (int n = 0; n < 4; n++) infile >> dummy >> dummy;
	infile >> dummy >> fields;
	infile >> dummy >> centerField;
	infile >> dummy >> fullFieldxBins;
	infile >> dummy >> fullFieldyBins;
	infile >> dummy >> downSamplingFactor;
	infile >> dummy >> approximatelyDownsampleNoiseCovariance;
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

	if (approximatelyDownsampleNoiseCovariance){
		sFullField->xBins = fullFieldxBins;
		sFullField->yBins = fullFieldyBins;
	} else {
		sFullField->xBins = xBins;
		sFullField->yBins = yBins;
	}

	sFullField->fBins = fBins;
	sFullField->fLength = fLength;
	sFullField->fStart = fStart;

	//Load in noise specs
	infile.open(NSpecsFilename.c_str(),fstream::in);
	infile >> dummy >> systemTemperatureConstant;
	infile >> dummy >> systemTemperatureCoeff;
	infile >> dummy >> systemTemperautreExp;
	infile >> dummy >> nAntennas;
	infile >> dummy >> observationTime;
	infile >> dummy >> useEmpiricalSigma;
	infile >> dummy >> baseRMSonInnerPixels;
	infile >> dummy >> throwOutModesWithObservationTimeRationsGreaterThanThis;
	infile >> dummy >> kernelConvolutionLimit;
	infile >> dummy >> throwOutModesThisMuchWorseThanTheBest;

	infile.close();	
}

//Loads a file in my format into a CVector
CVector loadFileIntoCVector(int field, int slice, string type, bool fullField){
	stringstream ss;
	if (fullField){
		ss << bigCubeDirectory << type;
	} else {
		ss << cubeDirectory << type;
	}
	if (!fullField) ss << dataCubeString1 << field;
	ss << datacubeString2 << slice << datacubeString3;
	fstream infile((ss.str()).c_str(),fstream::in);
	if (!fullField){
		CVector loadedVector(s);
		double value;
		for (int n = 0; n < nElements; n++){
			infile >> value;
			loadedVector.real[n] = value;
		}
		return loadedVector;
	} else {
		CVector loadedVector(sFullField);
		double value;
		for (int n = 0; n < nElements*fields; n++){
			infile >> value;
			loadedVector.real[n] = value;
		}
		return loadedVector;
	}
}

//The average of the synthesized beams for the two time slices is taken to be the synthesized beam from which we determine relative uv coverage
CVector loadAverageUVweights(bool fullField){
	if (fullField && approximatelyDownsampleNoiseCovariance) {
		return (loadFileIntoCVector(centerField, 0, "psf",fullField));
	} else {
		CVector uvweights0 = loadFileIntoCVector(centerField, 0, "psf",false);
		CVector uvweights1 = loadFileIntoCVector(centerField, 1, "psf",false);
		CVector uvweightsAvg = uvweights0 + uvweights1;
		for (int n = 0; n < uvweightsAvg.nElements; n++) uvweightsAvg.real[n] /= 2;
		return uvweightsAvg;
	}
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
	vector< vector<double> > avgObsTimes(sFullField->xBins, vector<double>(sFullField->yBins,0));
	for (int i = 0; i < sFullField->xBins; i++){
		for (int j = 0; j < sFullField->yBins; j++){
			for (int k = 0; k < sFullField->fBins; k++){
				avgObsTimes[i][j] += observationTimes[i][j][k] / (sFullField->fBins - nMaskedChannels); //averaged observation times
			}
		}
	}
	
	ofstream outfile;
	string ObservationTimesFilename = "obsTimes.dat";
	outfile.open(ObservationTimesFilename.c_str(), ios::trunc);	
	for (int j = 0; j < sFullField->yBins; j++){
		for (int i = 0; i < sFullField->xBins; i++) outfile << avgObsTimes[i][j] << " ";
		outfile << endl;
	}
	outfile.close();		
}

vector< vector< vector<double> > > calculateObservationTimes(CVector& uvweights, int nMaskedChannels){

	//Convert to uvk basis (Fourier perpendicular, real parallel)
	CVector FTWeights = uvweights.ijk2uvk();
	for (int n = 0; n < uvweights.nElements; n++){
		FTWeights.real[n] = sqrt(pow(FTWeights.real[n],2) + pow(FTWeights.imag[n],2));
		FTWeights.imag[n] = 0.0;
	}
	vector< vector< vector<double> > > obsTimes(uvweights.xBins, vector< vector<double> >(uvweights.yBins,vector<double>(uvweights.fBins,0)));
	for (int i = 0; i < uvweights.xBins; i++){
		for (int j = 0; j < uvweights.yBins; j++){
			for (int k = 0; k < uvweights.fBins; k++){
				obsTimes[i][j][k] = FTWeights.real[i*uvweights.yBins*uvweights.fBins + j*uvweights.fBins + k];
			}
		}
	}

	vector< vector<double> > avgObsTimes(uvweights.xBins, vector<double>(uvweights.yBins,0));
	for (int i = 0; i < uvweights.xBins; i++){
		for (int j = 0; j < uvweights.yBins; j++){
			for (int k = 0; k < uvweights.fBins; k++){
				avgObsTimes[i][j] += FTWeights.real[i*uvweights.yBins*uvweights.fBins + j*uvweights.fBins + k] / (uvweights.fBins - nMaskedChannels); //averaged observation times
			}
		}
	}
	
	//Normalize and convert to seconds
	double obsTimesSum = 0;
	for (int i = 0; i < uvweights.xBins; i++){
		for (int j = 0; j < uvweights.yBins; j++) obsTimesSum += avgObsTimes[i][j];
	}
	double maxObsTime = 0;
	double converstionToSeconds = observationTime * 3600.0 * nAntennas * (nAntennas - 1) / obsTimesSum;
	for (int i = 0; i < uvweights.xBins; i++){
		for (int j = 0; j < uvweights.yBins; j++){
			for (int k = 0; k < uvweights.fBins; k++){
				obsTimes[i][j][k] *= converstionToSeconds;
				if (obsTimes[i][j][k] > maxObsTime) maxObsTime = obsTimes[i][j][k];
			}
			avgObsTimes[i][j] *= converstionToSeconds;
			if (ignoreTopAndLeftEdge){
				if (i == 0 || j == 0){
					for (int k = 0; k < uvweights.fBins; k++){
						obsTimes[i][j][k] = 0;
					}	
					avgObsTimes[i][j] = 0;
				}
			}
		}
	}

	// Set Origin ObsTimes to 0
	for (int k = 0; k < uvweights.fBins; k++){
		obsTimes[uvweights.xBins/2][uvweights.yBins/2][k] = 0;	
	}

	//Delete very small observation times;
	for (int i = 0; i < uvweights.xBins; i++){
		for (int j = 0; j < uvweights.yBins; j++){
			for (int k = 0; k < uvweights.fBins; k++){
				if (obsTimes[i][j][k] < 1e-6 * maxObsTime) obsTimes[i][j][k] = 0;
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

//For a given slice, loads in all temps
vector<CVector> loadAllFields(int slice){
	vector<CVector> allTemps;
	for (int f = 0; f < fields; f++){
		CVector tempsHere = loadFileIntoCVector(f, slice, "temps",false);
		allTemps.push_back(tempsHere);
	}
	return allTemps;
}





double calculateOmegaPrime(vector<double>& freqs){	
	return 0.3948 * pow((fStart + (freqs[fBins -1] - fStart)/2)/150,-1.454);
}

void updateN(int n, double deltaNu, double OmegaPix, double OmegaPrime, double f, vector< vector< vector<double> > >& observationTimes, vector< vector < vector<double> > >& noiseCovarianceMatrixDiagonal){
	double Tsys = systemTemperatureConstant + systemTemperatureCoeff * pow(c/1e6/f,systemTemperautreExp);
	//cout << "For slice " << n << ", Tsys = " << Tsys << endl;
	double prefactor = pow(Tsys,2) * OmegaPrime / 2 / (deltaNu*1e6) / OmegaPix; //2 in the denominator matches Parsons et al. 2012 for two polarizations.
	for (int i = 0; i < observationTimes.size(); i++){
		for (int j = 0; j < observationTimes[0].size(); j++){
			if (observationTimes[i][j][n] == 0 || i == 0 || j == 0) noiseCovarianceMatrixDiagonal[i][j][n] = -1;
			else noiseCovarianceMatrixDiagonal[i][j][n] = prefactor / observationTimes[i][j][n];
		}
	}
}

void convertToNaturalWeighting(vector< vector < vector<double> > >& noiseCovarianceMatrixDiagonalFullField, CVector& fullFieldAverageUVweights){
	CVector FTWeights = fullFieldAverageUVweights.ijk2uvk();
	for (int n = 0; n < FTWeights.nElements; n++){
		FTWeights.real[n] = sqrt(pow(FTWeights.real[n],2) + pow(FTWeights.imag[n],2));
		FTWeights.imag[n] = 0.0;
	}

	for (int k = 0; k < FTWeights.fBins; k++){
		double FTWeightsAvg = 0;
		for (int i = 0; i < FTWeights.xBins; i++){
			for (int j = 0; j < FTWeights.yBins; j++){
				FTWeightsAvg += FTWeights.real[k + j*FTWeights.fBins + i*FTWeights.fBins*FTWeights.yBins];
			}
		}
		FTWeightsAvg /= (FTWeights.xBins * FTWeights.yBins);
		for (int i = 0; i < FTWeights.xBins; i++){
			for (int j = 0; j < FTWeights.yBins; j++){
				if (noiseCovarianceMatrixDiagonalFullField[i][j][k] > 0){ //TODO: test if this is good enough to catch all problems.
					noiseCovarianceMatrixDiagonalFullField[i][j][k] *= FTWeights.real[k + j*FTWeights.fBins + i*FTWeights.fBins*FTWeights.yBins] / FTWeightsAvg;
				} else {
					noiseCovarianceMatrixDiagonalFullField[i][j][k] = 0;
				}
			}
		}
	}
	return;
}

vector< vector < vector<double> > > CropAndDownSampleCovariance(vector< vector < vector<double> > >& noiseCovarianceMatrixDiagonalFullField){
	//Calculate the convolution kernel
	CVector croppingMask(sFullField);
	for (int i = 0; i < sFullField->xBins; i++){
		for (int j = 0; j < sFullField->yBins; j++){
			for (int k = 0; k < sFullField->fBins; k++){
				if (i < xBins*downSamplingFactor && j < yBins*downSamplingFactor) croppingMask.real[i*(sFullField->yBins)*(sFullField->fBins) + j*(sFullField->fBins) + k] = 1;
			}
		}
	}
	CVector croppingMaskFT = croppingMask.ijk2uvk();

	cout << "Cropping mask made..." << endl;
	
	vector< vector<double> > kernel(sFullField->xBins, vector<double>(sFullField->yBins,0));
	double kernelSum = 0;
	double kernelMax = 0;
	for (int i = 0; i < sFullField->xBins; i++){
		for (int j = 0; j < sFullField->yBins; j++){
			int here = i*(sFullField->yBins)*(sFullField->fBins) + j*(sFullField->fBins);
			kernel[i][j] = pow(croppingMaskFT.real[here],2) + pow(croppingMaskFT.imag[here],2);
			kernelSum += kernel[i][j];
			if (kernel[i][j] > kernelMax) kernelMax = kernel[i][j];
		}
	}

	//Eliminate small entries in the kernel, convert into sparse matrix form, and renormalize so that sum(kernel) = 1;
	vector<int> xCoords;
	vector<int> yCoords;
	vector<double> sparseKernel;
	double n = 0;
	for (int i = 0; i < sFullField->xBins; i++){
		for (int j = 0; j < sFullField->yBins; j++){
			if (kernel[i][j] / kernelMax > kernelConvolutionLimit){
				xCoords.push_back(i);
				yCoords.push_back(j);
				sparseKernel.push_back(kernel[i][j]/kernelSum);
			}
		} 
	}

	cout << "Kernel made..." << endl;

	//Apply convolution kernel to downsample
	int linearFieldsFactor = int(sqrt(fields));
	vector< vector < vector<double> > > croppedDiagN(downSamplingFactor*xBins, vector< vector<double> >(downSamplingFactor*yBins, vector<double>(fBins,0)));
	for (int i = 0; i < downSamplingFactor*xBins; i++){
		for (int j = 0; j < downSamplingFactor*yBins; j++){
			for (int n = 0; n < sparseKernel.size(); n++){
				int xCoordHere = (xCoords[n] - (sFullField->xBins)/2 + i*linearFieldsFactor + sFullField->xBins) % (sFullField->xBins);
				int yCoordHere = (yCoords[n] - (sFullField->yBins)/2 + j*linearFieldsFactor + sFullField->yBins) % (sFullField->yBins);
				for (int k = 0; k < fBins; k++){
					croppedDiagN[i][j][k] += sparseKernel[n] * noiseCovarianceMatrixDiagonalFullField[xCoordHere][yCoordHere][k];
				}
			}
		}
	}

	cout << "Kernel applied..." << endl;


	vector< vector <vector<double> > > croppedAndDownSampledDiagN(xBins, vector< vector<double> >(yBins, vector<double>(fBins,0)));
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			for (int k = 0; k < fBins; k++) {
				croppedAndDownSampledDiagN[i][j][k] = croppedDiagN[i+xBins*downSamplingFactor/2-xBins/2][j+yBins*downSamplingFactor/2-yBins/2][k] / downSamplingFactor / downSamplingFactor;
			}
		}
	}

	cout << "Downsampling done..." << endl;
	
	return croppedAndDownSampledDiagN;
}

void convertBackToUniformWeighting(vector< vector < vector<double> > >& noiseCovarianceMatrixDiagonal, CVector& averageUVweights){
	CVector FTWeights = averageUVweights.ijk2uvk();
	for (int n = 0; n < nElements; n++){
		FTWeights.real[n] = sqrt(pow(FTWeights.real[n],2) + pow(FTWeights.imag[n],2));
		FTWeights.imag[n] = 0.0;
	}

	for (int k = 0; k < fBins; k++){
		double FTWeightsAvg = 0;
		for (int i = 0; i < xBins; i++){
			for (int j = 0; j < yBins; j++){
				FTWeightsAvg += FTWeights.real[k + j*fBins + i*fBins*yBins];
			}
		}
		FTWeightsAvg /= (xBins * yBins);
		for (int i = 0; i < FTWeights.xBins; i++){
			for (int j = 0; j < FTWeights.yBins; j++){
				if (FTWeightsAvg < 1e-10 || FTWeights.real[k + j*fBins + i*fBins*yBins] < 1e-10 || noiseCovarianceMatrixDiagonal[i][j][k] == 0){ //TODO: double check this 
					noiseCovarianceMatrixDiagonal[i][j][k] = -1;
				} else {
					noiseCovarianceMatrixDiagonal[i][j][k] *= FTWeightsAvg / FTWeights.real[k + j*fBins + i*fBins*yBins];
				}
			}
		}
	}
	return;
}

void maskInconsistentUVCells(vector< vector< vector<double> > >& noiseCovarianceMatrixDiagonal){
	//Now throw out modes where the eigenvalue ratio of N is very large
	int highRatioModesThrownOut = 0;
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			double smallestNEig = 1e100;
			double largestNEig = 0;
			for (int k = 0; k < fBins; k++){
				if (noiseCovarianceMatrixDiagonal[i][j][k] > 0){
					if (noiseCovarianceMatrixDiagonal[i][j][k] > largestNEig) largestNEig = noiseCovarianceMatrixDiagonal[i][j][k];
					if (noiseCovarianceMatrixDiagonal[i][j][k] < smallestNEig && noiseCovarianceMatrixDiagonal[i][j][k] > 0) smallestNEig = noiseCovarianceMatrixDiagonal[i][j][k];
				}
			}
			double dist = sqrt(1.0*pow((i-xBins/2),2) + 1.0*pow((j-yBins/2),2));
			if (largestNEig/smallestNEig > throwOutModesWithObservationTimeRationsGreaterThanThis /*&& dist > baseRMSonInnerPixels*/){
				//cout << largestNEig/smallestNEig << endl;
				highRatioModesThrownOut++;
				for (int k = 0; k < fBins; k++) noiseCovarianceMatrixDiagonal[i][j][k] = -1;
			}
		}
	}
	cout << highRatioModesThrownOut << " uv modes thrown out due to highly inconsistant observation times (ratio > " << throwOutModesWithObservationTimeRationsGreaterThanThis << ")." << endl << endl;
}

void throwOutLeastObservedModes(vector< vector< vector<double> > >& noiseCovarianceMatrixDiagonal){
	double bestMode = 1e100;
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			for (int k = 0; k < fBins; k++){
				if ((noiseCovarianceMatrixDiagonal[i][j][k] > 0) && (noiseCovarianceMatrixDiagonal[i][j][k] < bestMode)) bestMode = noiseCovarianceMatrixDiagonal[i][j][k];
			}
		}
	}
	cout << bestMode << endl;
	int LoSthrownOut = 0;
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			double dist = sqrt(1.0*pow((i-xBins/2),2) + 1.0*pow((j-yBins/2),2));
			for (int k = 0; k < fBins; k++){
				if (noiseCovarianceMatrixDiagonal[i][j][k] > bestMode*throwOutModesThisMuchWorseThanTheBest){
					for (int k2 = 0; k2 < fBins; k2++) noiseCovarianceMatrixDiagonal[i][j][k2] = -1;
					LoSthrownOut++;
				}
				if (dist > baseRMSonInnerPixels) noiseCovarianceMatrixDiagonal[i][j][k] = -1;
			}
		}
	}
	cout << LoSthrownOut << " uv modes thrown out because at least one uvf cell is " << throwOutModesThisMuchWorseThanTheBest << " times worse than the best observed mode." << endl << endl;
}

void enforceCroppedFieldSymmetry(vector< vector< vector<double> > >& noiseCovarianceMatrixDiagonal){
	for (int k = 0; k < fBins; k++){
		for (int i = 0; i < xBins; i++){
			for (int j = 0; j < yBins; j++){
				if ((i == 0) || (j == 0) || ((i == xBins/2) && (j == yBins/2))){ 
					noiseCovarianceMatrixDiagonal[i][j][k] = -1;
				} else {
					double valueHere = noiseCovarianceMatrixDiagonal[i][j][k];
					double valueThere = noiseCovarianceMatrixDiagonal[xBins - i][yBins - j][k];
					if (valueHere == -1 || valueThere == -1){
						noiseCovarianceMatrixDiagonal[i][j][k] = -1;
						noiseCovarianceMatrixDiagonal[xBins - i][yBins - j][k] = -1;	
					} else {
						noiseCovarianceMatrixDiagonal[i][j][k] = (valueHere + valueThere)/2;	
						noiseCovarianceMatrixDiagonal[xBins - i][yBins - j][k] = (valueHere + valueThere)/2;	
					}
				}
			}
		}
	}
	return;
}

double calculateEmpiricalSigma(vector< vector < vector<double> > >& noiseCovarianceMatrixDiagonal){
	vector<CVector> Temps0 = loadAllFields(0);
	vector<CVector> Temps1 = loadAllFields(1);
	
	for (int f = 0; f < fields; f++){
		Temps0[f] = Temps0[f].ijk2uvk();
		Temps1[f] = Temps1[f].ijk2uvk();
	}

	int countOfUnmaskedBins = 0;
	for (int f = 0; f < fields; f++){
		for (int k = 0; k < fBins; k++){
			for (int i = 0; i < xBins; i++){
				for (int j = 0; j < yBins; j++){
					double dist = sqrt(1.0*pow((i-xBins/2),2) + 1.0*pow((j-yBins/2),2));
					if (noiseCovarianceMatrixDiagonal[i][j][k] == -1 || dist > baseRMSonInnerPixels){
						Temps0[f].real[i*yBins*fBins + j*fBins + k] = 0;
						Temps0[f].imag[i*yBins*fBins + j*fBins + k] = 0;
						Temps1[f].real[i*yBins*fBins + j*fBins + k] = 0;
						Temps1[f].imag[i*yBins*fBins + j*fBins + k] = 0;
					} else {
						countOfUnmaskedBins++;
					}
				}
			}
		}
	}
	for (int f = 0; f < fields; f++){
		Temps0[f] = Temps0[f].uvk2ijk();
		Temps1[f] = Temps1[f].uvk2ijk();
	}

	double totalEmpiricalSigmaSquared = 0;
	//cout << endl << "We have filtered out the top " << 100*fractionOfModesToThrowOut << "% noisiest modes and declared them 'unobserved'" << endl;
	for (int f = 0; f < fields; f++){
		double empricalSigmaSquared = 0;
		for (int k = 0; k < fBins; k++){
			double SigmaByFreq = 0;
			double theoreticalSigmaByFreq = 0;
			for (int i = 0; i < xBins; i++){
				for (int j = 0; j < yBins; j++){
					double temp0 = Temps0[f].real[i*yBins*fBins + j*fBins + k];
					double temp1 = Temps1[f].real[i*yBins*fBins + j*fBins + k];
					empricalSigmaSquared += pow(temp0-temp1,2);
					SigmaByFreq += pow(temp0-temp1,2);
					if (noiseCovarianceMatrixDiagonal[i][j][k] != -1) theoreticalSigmaByFreq += noiseCovarianceMatrixDiagonal[i][j][k];
				}
			}
			//cout << "Field " << f << ", Channel " << k << " Theoretical : " << sqrt(theoreticalSigmaByFreq / (xBins * yBins)) << endl;
			//cout << "Field " << f << ", Channel " << k << ": " << sqrt(SigmaByFreq / (xBins * yBins)) << endl; 
		}
		empricalSigmaSquared /= (countOfUnmaskedBins/fields);

		empricalSigmaSquared /= 2; //convert from RMS of the difference to sigma^2
		cout << "Emprical sigma for Field " << f << " is " << sqrt(empricalSigmaSquared) << endl;
		totalEmpiricalSigmaSquared += empricalSigmaSquared / fields;
	}
	cout << "Overall sigma: " << sqrt(totalEmpiricalSigmaSquared) << endl << endl;

	return sqrt(totalEmpiricalSigmaSquared);
}

void enforceSigma(vector< vector < vector<double> > >& noiseCovarianceMatrixDiagonal, vector< vector< vector<double> > >& observationTimes){
	double oldSigmaSquared = 0;
	int maskedElementsCounter = 0;
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			for (int k = 0; k < fBins; k++){
				double dist = sqrt(1.0*pow((i-xBins/2),2) + 1.0*pow((j-yBins/2),2));
				if (noiseCovarianceMatrixDiagonal[i][j][k] > 0 /*&& dist <= baseRMSonInnerPixels*/) oldSigmaSquared += noiseCovarianceMatrixDiagonal[i][j][k];
				else maskedElementsCounter++;
			}
		}
	}

	oldSigmaSquared /= (xBins*yBins*fBins - maskedElementsCounter);
	cout << "Theoretical Sigma: " << sqrt(oldSigmaSquared) << endl;
	if (useEmpiricalSigma){
		double empiricalSigma = calculateEmpiricalSigma(noiseCovarianceMatrixDiagonal);
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
	//Load in parameters and data
	cout << endl << "Now calculating the noise covariance matrix..." << endl;
	s = new Specs();
	sFullField = new Specs();
	loadSpecs("../../Specifications/dataSpecs.txt","../../cubeParameters.txt","../../Specifications/NSpecs.txt");
	cout << downSamplingFactor << endl;
	CVector fullFieldAverageUVweights = loadAverageUVweights(true);
	cout << "UVweights Loaded..." << endl;
	cout << xBins << yBins << fBins << endl;

	//Calculate observation times
	int nMaskedChannels = countMaskedChannels(fullFieldDataCubeMaskFilename);
	vector< vector< vector<double> > > observationTimes = calculateObservationTimes(fullFieldAverageUVweights,nMaskedChannels);
	printObsTimesToFile(observationTimes,nMaskedChannels);
	cout << "Observation Times Calculated..." << endl;

	//Calculate wavelengths and geometric factors
	vector<double> freqs = listFrequencies();
	double deltaF = freqs[0] - freqs[1];
	double bandwidth = fBins * deltaF; 
	double omegaPix = calculateOmegaPix(freqs);
	double OmegaPrime = calculateOmegaPrime(freqs); //Note: the f depenedence of omegaPix and OmegaPrime are assumed to cancel
	cout << "Omega' = " << OmegaPrime << " assuming a power law fit to frequency at " << (fStart + (freqs[fBins -1] - fStart)/2) << " MHz." << endl;

	//Calculate N
	vector< vector < vector<double> > > noiseCovarianceMatrixDiagonalFullField(sFullField->xBins, vector< vector<double> >(sFullField->yBins, vector<double>(sFullField->fBins,0)));
	for (int n = 0; n < fBins; n++) updateN(n,deltaF,omegaPix,OmegaPrime,freqs[n],observationTimes,noiseCovarianceMatrixDiagonalFullField);
	cout << "Full Field N calculated..." << endl;

	vector< vector < vector<double> > > noiseCovarianceMatrixDiagonal(s->xBins, vector< vector<double> >(s->yBins, vector<double>(s->fBins,0)));
	if (approximatelyDownsampleNoiseCovariance){
		//Convert to natural weighting
		convertToNaturalWeighting(noiseCovarianceMatrixDiagonalFullField, fullFieldAverageUVweights);
		cout << "Converted to Natural Weighting..." << endl;

		//Downsample and Crop N in real space
		noiseCovarianceMatrixDiagonal = CropAndDownSampleCovariance(noiseCovarianceMatrixDiagonalFullField);
		
		//Convert back to uniform weighting
		CVector averageUVweights = loadAverageUVweights(false);
		convertBackToUniformWeighting(noiseCovarianceMatrixDiagonal, averageUVweights);
		cout << "Converted back to uniform weighting..." << endl;
	} else {
		noiseCovarianceMatrixDiagonal = noiseCovarianceMatrixDiagonalFullField;
	}

	//Redunce computational complexity of the algorithm and eliminate imaginary components in real space
	maskInconsistentUVCells(noiseCovarianceMatrixDiagonal);
	throwOutLeastObservedModes(noiseCovarianceMatrixDiagonal);
	enforceCroppedFieldSymmetry(noiseCovarianceMatrixDiagonal);
	cout << "Masked appropriately..." << endl;

	//Rescale N to match the data	
	enforceSigma(noiseCovarianceMatrixDiagonal, observationTimes);
	printNtoFile(noiseCovarianceMatrixDiagonal);


	cout << "Noise covariance calculation complete." << endl << endl;
	return 0;
}


