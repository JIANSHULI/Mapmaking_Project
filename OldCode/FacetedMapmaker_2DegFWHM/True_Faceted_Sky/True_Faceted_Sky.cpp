#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <sstream>
#include "healpix_map.h" 
#include "healpix_data_io.h" 
#include "gsl/gsl_fit.h"
#include "healpix_map_fitsio.h" 
#include <xcomplex.h>
#include "vec3.h"
#include "pointing.h"
#include "lsconstants.h"
#include "trafos.h"
#include <stdio.h>
#include <string.h>

using namespace std;

//Global variables
const double c = 299792458;
double south, east, up, freq, startingLST, LST, duration, timeStep, latitude, longitude, xpolOrientationDegreesEastofNorth, UTCtime, facetRA, facetDec, angularResolution, facetSize;
int res, PSFextensionBeyondFacetFactor;
string specsFile = "./True_Faceted_Sky/Specifications.txt";
string facetSpecsFile, dataProductFilePrefix;
string componentsFile = "components.dat";
int componentsToFit = 9;
string FITS_directory, alm_directory, GSMformat1, GSMformat2, polarization, dateString, timeString, resultsFolder;
int startingPower, endingPower;
int principalComps = 3; //Hardcoded...won't work otherwise
bool beamIsOne = false;
int nAlphaCropped, nDeltaCropped;
double arrayLat, latInRad, facetDecInRad, facetRAinRad, angResInRad, k;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// USEFUL DATA STRUCTURES
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Stores complex numbers (e.g. visibilities). Can multiply them by real (double precision) or other complex numbers
struct complex{
	double re, im;
	complex(){}
	complex(double reIn, double imIn) : re(reIn), im(imIn) {}
	complex conj(){ return complex(re, -im); }
	complex operator*(const complex& that){ return complex(re*that.re - im*that.im, re*that.im + im*that.re); }
	complex operator*(const double that){ return complex(re*that, im*that); }
	complex operator+(const complex& that){ return complex(re + that.re, im + that.im); }
	complex operator-(const complex& that){ return complex(re - that.re, im - that.im); }
	void print(){ cout << re << " + " << im << "i"; }
};


struct horizPoint; //Predeclaration

//When appropriate, +x is east, +y is north,  +z is up
struct cartVec{
	double x,y,z;	
	cartVec(){}
	cartVec(double xIn, double yIn, double zIn) : x(xIn), y(yIn), z(zIn) {}
	cartVec operator+(const cartVec& that){ return cartVec(x+that.x, y+that.y, z+that.z); }
	cartVec operator-(const cartVec& that){ return cartVec(x-that.x, y-that.y, z-that.z); }
	cartVec operator*(const double that){ return cartVec(x*that, y*that, z*that); }
	double dot(cartVec that){
		return (x*that.x + y*that.y + z*that.z);
	}
	cartVec cross(cartVec that){
		return cartVec(y*that.z - z*that.y, z*that.x - x*that.z, x*that.y - y*that.x);
	}
	cartVec normalize(){
		return (cartVec(x,y,z) * (1.0 / sqrt(x*x + y*y + z*z)));
	}
	horizPoint toHoriz();
	void print(){ cout << "["<< x << ", " << y << ", " << z << "]"; }
};

//Alt is 0 at the horizon, pi/2 at the zenith. Az is radians east of north
struct horizPoint{
	double alt, az;
	horizPoint(){}
	horizPoint(double altIn, double azIn) : alt(altIn), az(azIn) {}
	cartVec toVec(){
		cartVec cCoord(sin(az)*cos(alt), cos(az)*cos(alt), sin(alt));
		return cCoord;
	}
	double greatCircleAngle(horizPoint hp){
		horizPoint thisHP(alt, az);
		cartVec here = thisHP.toVec();
		cartVec there = hp.toVec();
		double dist = sqrt(pow(here.x - there.x,2) + pow(here.y - there.y,2) + pow(here.z - there.z,2));
		return 2*asin(dist/2);
	}
};

//Converts cartesian vectors to altitude and azimuth. I needed to delcate it outside the struct because horizPoint hadn't been declared yet.
horizPoint cartVec::toHoriz(){
	return horizPoint(asin(z/sqrt(x*x + y*y + z*z)), fmod(atan2(x,y)+2*pi,2*pi)); 
}

//Creates an equatorial pointing (RA and Dec) which can be converted into a horizontal pointing given an LST (assuming an array latitude as a global variable)
struct equaPoint{
	double ra, dec;
	equaPoint(){}
	equaPoint(double raIn, double decIn) : ra(raIn), dec(decIn) {}
	horizPoint toHoriz(double LST){
		double lha = pi/12.0*LST - ra; //radians
		horizPoint hp;
    	hp.alt = asin(sin(latInRad) * sin(dec) + cos(latInRad) * cos(dec) * cos(lha));
    	hp.az = atan2(sin(lha) * cos(dec), cos(lha) * cos(dec) * sin(latInRad) - sin(dec) * cos(latInRad)) + pi;
    	return hp;
	}
};



//Loads the relevant information about the GSM files and the observation into global variables
void loadSpecs(){
	fstream infile(specsFile.c_str(),fstream::in);
	string dummy;
	infile >> dummy >> FITS_directory;
	infile >> dummy >> GSMformat1;
	infile >> dummy >> GSMformat2;
	infile >> dummy >> res;
	infile.close();

	fstream infile2(facetSpecsFile.c_str(),fstream::in);
	infile2 >> dummy >> dummy >> freq;
	infile2 >> dummy >> dummy >> dummy;
	infile2 >> dummy >> dummy >> dummy;
	infile2 >> dummy >> dummy >> dummy;
	infile2 >> dummy >> dummy >> dummy;
	infile2 >> dummy >> dummy >> dummy;
	infile2 >> dummy >> dummy >> dummy;
	infile2 >> dummy >> dummy >> dummy;
	infile2 >> dummy >> dummy >> dummy;
	infile2 >> dummy >> dummy >> dataProductFilePrefix;
	infile2 >> dummy >> dummy >> dummy;
	infile2 >> dummy >> dummy >> dummy;
	infile2 >> dummy >> dummy >> facetRA;
	infile2 >> dummy >> dummy >> facetDec;
	infile2 >> dummy >> dummy >> facetSize;
	infile2 >> dummy >> dummy >> angularResolution;
	infile2 >> dummy >> dummy >> PSFextensionBeyondFacetFactor;
	infile2.close();
	
	nAlphaCropped = int(round(facetSize / angularResolution * PSFextensionBeyondFacetFactor / 2)*2);
	nDeltaCropped = int(round(facetSize / angularResolution * PSFextensionBeyondFacetFactor / 2)*2);
	cout << nAlphaCropped << " " << nDeltaCropped << endl;
	arrayLat = facetDec;
	facetDecInRad = facetDec*pi/180;
	facetRAinRad = facetRA*pi/180;
	angResInRad = angularResolution*pi/180;
	latInRad = pi/180.0*arrayLat;
	LST = 24.0/360*facetRA;
}


//Loads in all 3 principal components of the GSM and return them as a vector.
vector< Healpix_Map<double> > loadGSMpcomps(){
	vector< Healpix_Map<double> > GSMpcomps;
	for (int p = 1; p <= principalComps; p++){
		stringstream ss;
		ss << FITS_directory << GSMformat1 << p << GSMformat2 << res;
		string GSMfile = ss.str();
		cout << "Now loading " << GSMfile << endl;
		Healpix_Map<double> map; 
		read_Healpix_map_from_fits(GSMfile,map); 
		GSMpcomps.push_back(map);
	}
	return GSMpcomps;
}

//Load the components data and properly reweight and rescale the GSM principal componenets
Healpix_Map<double> computeGSM(vector< Healpix_Map<double> >& GSMpcomps){
	//Load in the componenets data
	stringstream ss;
	ss << FITS_directory << componentsFile;
	fstream infile((ss.str()).c_str(),fstream::in);
	vector<double> pcompFreqs(componentsToFit,0);
	vector<double> pcompTemps(componentsToFit,0);
	vector< vector<double> > weights(3, vector<double>(componentsToFit,0));
	for (int n = 0; n < componentsToFit; n++){
		double f,T,p1,p2,p3;
		infile >> f >> T >> p1 >> p2 >> p3;
		pcompFreqs[n] = f;
		pcompTemps[n] = T;
		weights[0][n] = p1;
		weights[1][n] = p2;
		weights[2][n] = p3;
	}
	infile.close();
	
	//Determine the temperature though a power law fit
	double* logFreqs = new double[componentsToFit];
	double* logTemps = new double[componentsToFit];
	for (int n = 0; n < componentsToFit; n++){
		logFreqs[n] = log10(pcompFreqs[n]);
		logTemps[n] = log10(pcompTemps[n]);
	}
	double c0, c1, cov00, cov01, cov11, sumsq;
	gsl_fit_linear(logFreqs, 1, logTemps, 1, componentsToFit, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
	delete[] logFreqs;
	delete[] logTemps;

	double overallTemperature = pow(10.0,c0)*pow(freq,c1);
	cout << "A f^" << c1 << " power law yields a GSM temperature normalization: " << overallTemperature << " K." << endl;

	//Determine the weights through linear interpolation
	int freqCounter = 1;
	while (freqCounter < componentsToFit){
		if (pcompFreqs[freqCounter] > freq) break;
		freqCounter++;
	}
	vector<double> w(3,0);
	for (int i = 0; i < 3; i++){
		double slope = (weights[i][freqCounter] - weights[i][freqCounter-1]) / (pcompFreqs[freqCounter] - pcompFreqs[freqCounter-1]);
		w[i] = weights[i][freqCounter-1] + slope*(freq - pcompFreqs[freqCounter-1]);
	} 

	//Combine the components and return the final map;
	Healpix_Map<double> GSM = GSMpcomps[0];
	for (int n = 0; n < 12*res*res; n++) GSM[n] = (GSMpcomps[0][n]*w[0] + GSMpcomps[1][n]*w[1] + GSMpcomps[2][n]*w[2]) * overallTemperature;
	return GSM;
}


//Converts vector direction of x (toward Gal Center), y, z (toward NGP) to (altitude, azimuth east from north) at the specified location and local siderial time
equaPoint convertGalToEquaPoint(pointing& galPointing){
	double b = pi/2.0 - galPointing.theta;
	double l = galPointing.phi;
	double pole_ra = 2.0*pi/360.0*192.859508;
    double pole_dec = 2.0*pi/360.0*27.128336;
    double posangle = 2.0*pi/360.0*(122.932-90.0);
    double raHere = atan2((cos(b)*cos(l-posangle)), (sin(b)*cos(pole_dec) - cos(b)*sin(pole_dec)*sin(l-posangle))) + pole_ra;
    double decHere = asin(cos(b)*cos(pole_dec)*sin(l-posangle) + sin(b)*sin(pole_dec));
    return equaPoint(raHere, decHere);
}


/*//This function figures out the unit vectors for the direction of the pixelization either in the plane of the array for the FFT or in the facet for the purposes of the primary beam calculation
//They are stored as principalAxesOfTheFacet[+RA, +Dec, +RA cross +Dec]
vector<cartVec> findFacetPrincipalAxes(){
	cartVec center = (equaPoint(facetRAinRad, facetDecInRad)).toHoriz(LST).toVec();
	cartVec plusDec = (equaPoint(facetRAinRad, facetDecInRad + 1e-6)).toHoriz(LST).toVec();
	cartVec plusRA = (equaPoint(facetRAinRad + 1e-6, facetDecInRad)).toHoriz(LST).toVec();

	cartVec	plusRADirection = plusRA - center;
	cartVec plusDecDirection = plusDec - center;
	cartVec RAcrossDecDirection = plusRADirection.cross(plusDecDirection);

	vector<cartVec> principalAxesOfTheFacet;
	principalAxesOfTheFacet.push_back(plusRADirection.normalize());
	principalAxesOfTheFacet.push_back(plusDecDirection.normalize());
	principalAxesOfTheFacet.push_back(RAcrossDecDirection.normalize());
	return principalAxesOfTheFacet;
}*/

//This is the main function that computes the map wrapped around the sphere
vector< vector<double> > computeFacetMap(Healpix_Map<double>& GSM){
	/*vector<cartVec> principalAxesOfTheFacet = findFacetPrincipalAxes();
	principalAxesOfTheFacet[0].print(); cout << endl;
	principalAxesOfTheFacet[1].print(); cout << endl;
	principalAxesOfTheFacet[2].print(); cout << endl;*/

	vector< vector<double> > facetMap(nAlphaCropped, vector<double>(nDeltaCropped,0.0));
	vector< vector<int> > facetCounts(nAlphaCropped, vector<int>(nDeltaCropped,0));
	int N = 12*res*res;
	for (int n = 0; n < N; n++){
		pointing thisGalPoint = GSM.pix2ang(n);
		equaPoint thisPoint = convertGalToEquaPoint(thisGalPoint);
		cartVec	thisVec = (thisPoint.toHoriz(LST)).toVec();
		if (thisVec.z > 0){ 
			cartVec thisVecExtended = thisVec * (1.0/thisVec.z);
			int alphaIndex = int(round(thisVecExtended.x / angResInRad));
			int deltaIndex = int(round(thisVecExtended.y / angResInRad));
			if ((alphaIndex >= -nAlphaCropped/2) && (alphaIndex < nAlphaCropped/2) && (deltaIndex >= -nDeltaCropped/2) && (deltaIndex < nDeltaCropped/2)){
				facetMap[alphaIndex + nAlphaCropped/2][deltaIndex + nDeltaCropped/2] = facetMap[alphaIndex + nAlphaCropped/2][deltaIndex + nDeltaCropped/2] + GSM[n];
				facetCounts[alphaIndex + nAlphaCropped/2][deltaIndex + nDeltaCropped/2] = facetCounts[alphaIndex + nAlphaCropped/2][deltaIndex + nDeltaCropped/2] + 1;
			}
		}
	}
	cout << endl << endl;
	for (int a = 0; a < nAlphaCropped; a++){
		for (int d = 0; d < nDeltaCropped; d++){
			if (facetCounts[a][d] > 0) facetMap[a][d] = facetMap[a][d]/facetCounts[a][d];
			cout << facetCounts[a][d] << " ";
		}
		cout << endl;
	}

	return facetMap;
}

//Saves the results to file with a very detailed filename
void saveResultsToFile(vector< vector<double> >& facetMap){
	ofstream outFile;
	outFile.precision(16);
	stringstream ss;
	ss << dataProductFilePrefix << "true_sky.dat";
	string outFileName = ss.str();
	cout << "Now writing to " << outFileName << endl;

	outFile.open(outFileName.c_str(),ios::trunc);
	for (int a = 0; a < nAlphaCropped; a++){
		for (int d = 0; d < nDeltaCropped; d++){
			outFile << facetMap[a][d] << " ";
		}
		outFile << endl;
	}
	outFile.close();
}

int main(int argc, char *argv[]){ 

	// Load in all parameters
	facetSpecsFile = argv[1];
	loadSpecs();
	vector< Healpix_Map<double> > GSMpcomps = loadGSMpcomps();
	
	//Loop over all frequencies
	k = 2.0 * pi * freq*1e6 / c;
	Healpix_Map<double> GSM = computeGSM(GSMpcomps);
	vector< vector<double> > facetMap = computeFacetMap(GSM);
	saveResultsToFile(facetMap);

	cout << "Program executed successfully." << endl << endl;
	return 0;
}

