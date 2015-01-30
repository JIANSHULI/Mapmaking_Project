/* 
TODO:
* Reevaluate the way the primary beam is projected onto the facet.
* Make sure I figure out the signs on the k.b terms.
* Look again at whether we can easily fit polynomials to the PSF as a function of position


DONE:
* Create a "true" map by projecting the GSM onto a tangent plane at the facet location, so I can see how a convolution with the PSF compares to my dirty map.
* Test to see how the number of integrations in converges on the true map.
* Examine HERA, Omniscope, and MWA instantaneous beams at the zenith using Max's formula
* Make sure my explicitly calculated A^t (the matrix that converts from baseline space to image space) gives the same result in both parts of the code (the part that makes a map using an FFT and the part that computes the PSF).
	This works as long as I compare facets where the e^-k.b has been done on a per-basline basis and not on a gridded baseline basis...which gives the wrong answer to a few percent.
		This discrepancy doesn't matter in the long run, since the gridded phase-shifting cancels out when I calculate A^t N^-1 A
* It might be possible to make the code run faster by only computing Fourier modes where there are observations in any given facet.
	* Yes, I can save a bunch of time in the most time-intensive step by not summing over modes with Ninv[i][j] = 0;
* I don't need to know how pixels outside the facet are calculated...perhaps I should only look at how pixels inside the facet are affected by pixels outside the facet.
*/


#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <string>
#include <cstring>
#include <sstream>
#include "fftw3.h"
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <sstream>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SETTINGS AND GLOBAL VARIABLES
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

string specsFilename = "Specifications.txt";

//Settings to load
double freq; //MHz
string polarization;
string baselinesFile, visibilitiesFolder, unnormalizedFullMapFilename, finalMapFilename, PSFfilename, noiseCovFilename, DmatrixFilename, dataProductFilePrefix;
double arrayLat; //Degrees
double arrayLong; //Degrees
double facetRA; //Degrees
double facetDec; //Degrees;
double integrationTime; //seconds
double channelBandwidth; //MHz
double snapshotTime; //seconds
double angularResolution; //degrees
double facetSize; //degrees
int PSFextensionBeyondFacetFactor; //Should be an odd integer. 1 is no padding.
double maximumAllowedAngleFromPBCenterToFacetCenter = 10; //degrees
double xpolOrientationDegreesEastofNorth;
double noiseStd; //Jy on each visibilility...this number is totally made up
bool gaussianPrimaryBeam; //if true, use guassian PB at zenith, else use MWA dipoles (ala omniscope)
double primaryBeamFWHM; //FWHM in degrees if using a guassian PB

//Hard-coded settings and global variables
bool PSFpeaksAtOne = false; //Otherwise, the diagonal of the PSF matrix is set to be all ones
const double pi = 3.1415926535897932384626433832795;
const double c = 299792000; //m/s
int nAltBeamPoints = 45; //running from 0 to 90 in increments of 90/44
int nAzBeamPoints = 180; //running for 0 to 360 in increments of 360/179
double firstBeamFreq = 110; //MHz
double beamDeltaFreq = 10; //MHz
double lastBeamFreq = 190; //MHz
string beamDirectory = "../MWA_Primary_Beams/mwa_beam_";
int nFreqFiles = int((lastBeamFreq - firstBeamFreq)/10 + 1);
double deltaBeamAlt = 90.0 / (nAltBeamPoints-1);
double deltaBeamAz = 360.0 / (nAzBeamPoints-1);

//Global variables to compute later
int nIntegrationsToUse = 0;
int nAlphaFacet, nDeltaFacet, nAlphaCropped, nDeltaCropped, nBaselines, nLSTs, nSnapshots;
double k, facetDecInRad, facetRAinRad, angResInRad, latInRad, temperatureConversionFactor;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// USEFUL DATA STRUCTURES
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Stores 2D coordinates a and d
struct coord{
	int a, d;
	coord(int aIn, int dIn) : a(aIn), d(dIn) {}
};

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// LOADING AND INTIALIZATION FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Loads specifications from Specifications.txt or a specified filename;
void loadSpecs(){
	fstream infile(specsFilename.c_str(),fstream::in);
	string dummy;
	infile >> dummy >> dummy >> freq;
	infile >> dummy >> dummy >> polarization;
	infile >> dummy >> dummy >> baselinesFile;
	infile >> dummy >> dummy >> visibilitiesFolder;
	infile >> dummy >> dummy >> unnormalizedFullMapFilename;
	infile >> dummy >> dummy >> finalMapFilename;
	infile >> dummy >> dummy >> PSFfilename;
	infile >> dummy >> dummy >> noiseCovFilename;
	infile >> dummy >> dummy >> DmatrixFilename;
	infile >> dummy >> dummy >> dataProductFilePrefix;
	infile >> dummy >> dummy >> arrayLat;
	infile >> dummy >> dummy >> arrayLong;
	infile >> dummy >> dummy >> facetRA;
	infile >> dummy >> dummy >> facetDec;
	infile >> dummy >> dummy >> facetSize;
	infile >> dummy >> dummy >> angularResolution;
	infile >> dummy >> dummy >> PSFextensionBeyondFacetFactor;
	infile >> dummy >> dummy >> snapshotTime;
	infile >> dummy >> dummy >> integrationTime;
	infile >> dummy >> dummy >> channelBandwidth;
	infile >> dummy >> dummy >> gaussianPrimaryBeam;
	infile >> dummy >> dummy >> primaryBeamFWHM;
	infile >> dummy >> dummy >> xpolOrientationDegreesEastofNorth;
	infile >> dummy >> dummy >> maximumAllowedAngleFromPBCenterToFacetCenter;
	infile >> dummy >> dummy >> noiseStd;
	infile.close();

	//Other global variables to set
	nAlphaFacet = int(round(facetSize / angularResolution / 2)*2);
	nDeltaFacet = int(round(facetSize / angularResolution / 2)*2);
	nAlphaCropped = nAlphaFacet * PSFextensionBeyondFacetFactor;
	nDeltaCropped = nDeltaFacet * PSFextensionBeyondFacetFactor;
	k = 2 * pi * freq * 1e6 / c; 
	facetDecInRad = facetDec*pi/180;
	facetRAinRad = facetRA*pi/180;
	angResInRad = angularResolution*pi/180;
	latInRad = pi/180.0*arrayLat;
	temperatureConversionFactor = 2 * 1.3806488e-23 * pow(freq*1e6,2) / (c*c) * 1e26; //convert Kelvin to Jy
	unnormalizedFullMapFilename = dataProductFilePrefix + unnormalizedFullMapFilename;
	finalMapFilename = dataProductFilePrefix + finalMapFilename;
	PSFfilename = dataProductFilePrefix + PSFfilename;
	noiseCovFilename = dataProductFilePrefix + noiseCovFilename;
	DmatrixFilename = dataProductFilePrefix + DmatrixFilename;
}

//Loads baselines (only one orientation, no complex conjugates, from the file under the global variable "baselinesFile", which is in "south east up" format)
vector<cartVec> loadBaselines(vector<int>& baselineRedundancy){
	cout << "Now loading all the baseline vectors." << endl << "[NOTE: THIS ASSUMES THAT THE DATA FORMAT IS SOUTH EAST UP, THOUGH THE INTERNAL FORMAT IS EAST NORTH UP (RIGHT HANDED SYSTEM)]" << endl;
	vector<cartVec> baselines;
	double south, east, up;
	int multiplicity;
	fstream infile(baselinesFile.c_str(),fstream::in);
	cout << baselinesFile << endl;
	while(infile >> south >> east >> up >> multiplicity){
		cartVec thisBaseline;
		thisBaseline.y = -south; 
		thisBaseline.x = east;
		thisBaseline.z = up;
		baselines.push_back(thisBaseline);
		baselineRedundancy.push_back(multiplicity);
	}
	infile.close();
	nBaselines = baselines.size();
	return baselines;
}

//Loads the pre-computed principal axes of the plane of the array, which is, in general, slightly sloped relative to the the EW/NS plane, and a third axis perpendicular to the two and pointed mostly up
//The two axes must be perpendicular to one another. By convention, the first vector will be mostly East with no N/S component and the second will be mostly South
//The rotation can be paraterized by two cartesian rotations, one about the x (EW) axis and one about the y axis (SN)...this is completely general.
//Since rotations are not commutative, I'll first perform the rotation around the x axis and then the rotation about the y axis. This will preserve the property that x has no N/S component.
vector<cartVec> loadArrayPrincipalAxes(){
	cout << "Now loading array's principal axes...[NOTE: THESE ARE HARD CODED TO AN ARBITARY (SMALL) VALUE OR 0 FOR NOW.]" << endl;
	double thetaX = 0;//.02; //radians of rotation about the EW axis taking north toward up.
	double thetaY = 0;//.04; //radians of rotation about the NS axis taking east toward up.
	/*cartVec axisE(1, 0, 0); //Due east
	cartVec axisN(0, 1, 0); //Due south
	cartVec axisErotX(1, 0, 0);
	cartVec axisNrotX(0, cos(thetaX), sin(thetaX));*/
	cartVec axisErotXY(cos(thetaY), 0, sin(thetaY));
	cartVec axisNrotXY(-sin(thetaX)*sin(thetaY), cos(thetaX), sin(thetaX)*cos(thetaY)); 
	cartVec axisUrotXY = axisErotXY.cross(axisNrotXY);
	vector<cartVec> arrayPrincipalAxes;
	arrayPrincipalAxes.push_back(axisErotXY);		
	arrayPrincipalAxes.push_back(axisNrotXY);
	arrayPrincipalAxes.push_back(axisUrotXY);
	return arrayPrincipalAxes;
}

//This figures out the projection of the baselines into the best-fit plane of the array. This is needed for the 2D FFT.
//The projected baselines are written in terms of their components in terms of the array principal axes, not in terms of East, North, and Up.
vector<cartVec> calculateProjectedBaselines(vector<cartVec>& baselines, vector<cartVec>& arrayPrincipalAxes){
	vector<cartVec> projectedBaselines;
	for (int b = 0; b < nBaselines; b++){
		cartVec	projectedBaselineInENU = baselines[b] - arrayPrincipalAxes[2]*(baselines[b].dot(arrayPrincipalAxes[2]));
		cartVec projectedBaselineInXYZ(projectedBaselineInENU.dot(arrayPrincipalAxes[0]), projectedBaselineInENU.dot(arrayPrincipalAxes[1]), 0.0);
		projectedBaselines.push_back(projectedBaselineInXYZ);
	} 
	return	projectedBaselines;
}

//Loads visibilities, which have an LST and a real and imaginary component, from a file with the format I used for the omniscope
//Format is allVisibilites[baseline number same as baselines vector][different measurements].re/im
//It is assumed that all visibilities have the same set of LSTs, which is taken from the first baseline loaded
vector< vector<complex> > loadVisibilities(vector<cartVec>& baselines, vector<double>& LSTs){	
	vector< vector<complex> > allVisibilities;
	cout << "Now loading all visibilities for " << freq << " MHz and " << polarization << " polarization over " << nBaselines << " baselines..." << endl;
	for (int n = 0; n < nBaselines; n++){
		cout << " " << floor(100.0 * n / nBaselines) << "% done. \r" << std::flush;

		//Format filename
		stringstream ss;
		ss << visibilitiesFolder  << -baselines[n].y << "_m_south_" << baselines[n].x << "_m_east_" << baselines[n].z << "_m_up_" << polarization << "_pol_" << freq << "_MHz.dat";
		//ss << visibilitiesFolder << "Visibilties_for_" << round(-baselines[n].y) << "_m_south_" << round(baselines[n].x) << "_m_east_" << round(baselines[n].z) << "_m_up_" << polarization << "_pol_" << freq << "_MHz.dat";
		//ss << visibilitiesFolder << "Visibilties_for_" << -baselines[n].y << "_m_south_" << round(baselines[n].x) << "_m_east_" << round(baselines[n].z) << "_m_up_" << polarization << "_pol_" << freq << "_MHz.dat";
		string infilename = ss.str();
		cout << infilename << endl;
		fstream infile(infilename.c_str(),fstream::in);
		
		//Load all LST, re, and im into a vector of vectors of visibilities
		vector<complex> thisBaselinesVisibilities;
		double localSiderealTime, real, imaginary;
		while(infile >> localSiderealTime >> real >> imaginary){
			if (n == 0) LSTs.push_back(localSiderealTime);
			complex thisVisibility(real, imaginary);
			thisBaselinesVisibilities.push_back(thisVisibility);
		}
		infile.close();
		if (thisBaselinesVisibilities.size() == 0) cout << "WARNING: Cannot find " << infilename << endl;
		allVisibilities.push_back(thisBaselinesVisibilities);
	}
	nLSTs = LSTs.size();
	cout << "Done.                  " << endl;  
	return allVisibilities;
}

vector< vector<double> > loadNoiseVarianceOnEachVisibiltiy(vector<cartVec>& baselines, vector<double>& LSTs, vector<int> baselineRedundancy){
	//TODO: ENVENTUALLY I'LL WANT TO LOAD A MODEL FOR THE NOISE ON ANTENNA AND THEN COMPUTE WHAT THE PER VISIBILITY NOISE IS. FOR NOW I'LL JUST USE A SIMPLE MODEL.
	cout << "Now loading and computing the noise variance on each visibility...[NOTE: NOISE SET TO A SIMPLISTIC MODEL FOR NOW.]" << endl;
	vector< vector<double> > noiseVarianceOnEachVisibility(nBaselines, vector<double>(nLSTs,0));
	for (int t = 0; t < nLSTs; t++){
		for (int b = 0; b < nBaselines; b++){
			noiseVarianceOnEachVisibility[b][t] = pow(noiseStd,2)/(integrationTime * channelBandwidth * 1e6 * baselineRedundancy[b]);
		}
	}
	return noiseVarianceOnEachVisibility;
}

//This function loads the location of the primary beam in horizontal coordinates as a function of LST
vector<horizPoint> loadPBPointings(vector<double>& LSTs){
	//TODO: EVENTUALLY, WE WANT TO GENERALIZE THIS TO WHERE THE PRIMARY BEAM IS POINTED, BUT FOR NOW WE'LL ASSUME IT'S THE ZENITH
	cout << "Now loading all the primary beam pointings..." << endl;
	vector<horizPoint> PBpointings;
	for (int n = 0; n < nLSTs; n++){
		horizPoint thisPointing(pi/2, 0);
		PBpointings.push_back(thisPointing);
	}
	return PBpointings;
}

//This funciton loads the primary beam into a array of discretized values of alt and az.
//The data file, somewhat inconsistently, has azimuth 0 in the direction pointed by the XX polarization and increases CCW
vector< vector<double> > loadDiscretizedPrimaryBeam(){
	cout << "Now loading the primary beams for the nearest two frequencies and interpolating/extrapolating between them..." << endl;
	cout << "[NOTE: FOR NOW, PRIMARY BEAM IS ASSUMED TO BE THE SAME FOR ALL OBSERVATIONS]" << endl;
	vector< vector<double> > primaryBeamDiscretized(nAltBeamPoints, vector<double>(nAzBeamPoints, 0.0));
	double freq1 = -1;
	double freq2 = -1;
	for (double f = firstBeamFreq; f <= lastBeamFreq; f+=beamDeltaFreq){
		if (fabs(f - freq) < fabs(freq1 - freq)) freq1 = f;	
	} 
	for (double f = firstBeamFreq; f <= lastBeamFreq; f+=beamDeltaFreq){
		if ((fabs(f - freq) < fabs(freq2 - freq)) && (f != freq1)) freq2 = f;
	}
	stringstream freq1stream, freq2stream;
	freq1stream << beamDirectory << polarization << "_" << freq1 << ".dat";
	freq2stream << beamDirectory << polarization << "_" << freq2 << ".dat";
	string file1 = freq1stream.str();
	string file2 = freq2stream.str();
	fstream infile1(file1.c_str(),fstream::in);
	fstream infile2(file2.c_str(),fstream::in);
	double gain1, gain2;
	for (int alt = 0; alt < nAltBeamPoints; alt++){
		for (int az = 0; az < nAzBeamPoints; az++){
			infile1 >> gain1;
			infile2 >> gain2;
			primaryBeamDiscretized[alt][az] = (gain2 - gain1)/(freq2 - freq1) * (freq - freq1) + gain1;
		}
	}
	return primaryBeamDiscretized;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GEOMETRY FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//This function computes the altitude and azimuth of the facet ceneter for all LSTs observed.
vector<horizPoint> computeFacetCenterPointings(vector<double>& LSTs){
	cout << "Now computing the horizontal angle to the facet center for all LSTs..." << endl;
	vector<horizPoint> facetCenterPointings;
	equaPoint facetCenter(facetRAinRad, facetDecInRad);
	for (int n = 0; n < nLSTs; n++) facetCenterPointings.push_back(facetCenter.toHoriz(LSTs[n]));
	return facetCenterPointings;
}

//This function determines whether the distance between the facet center and the center of the maximumAllowedAngleFromPBCenterToFacetCenter
vector<bool> determineIntegrationsToUse(vector<horizPoint>& PBpointings, vector<horizPoint>& facetCenterPointings){
	vector<bool> toUse;
	for (int n = 0; n < PBpointings.size(); n++){ 
		toUse.push_back(facetCenterPointings[n].greatCircleAngle(PBpointings[n]) < pi/180*maximumAllowedAngleFromPBCenterToFacetCenter);
		if (toUse[n]) nIntegrationsToUse++;
	}
	cout << nIntegrationsToUse << " of the " << PBpointings.size() << " total integrations have the facet center within " << maximumAllowedAngleFromPBCenterToFacetCenter << " degrees of the primary beam center." << endl;
	return toUse;
}

//This functiona assigns LST indices to snapshot indices, which are gridded together
vector< vector<int> > assignLSTsToSnapshots(vector<bool> integrationsToUse){
	vector< vector<int> > snapshotLSTindices;
	for (int t = 0; t < nLSTs; t++){		
		if (integrationsToUse[t]){ //if ready to start a snapshot
			vector<int> integrationIndices;
			for (int t2 = 0; (t2 < snapshotTime/integrationTime && t < nLSTs); t2++){ //loop over a whole snapshot or until the end of observations
				if (integrationsToUse[t]) integrationIndices.push_back(t); //accumulate good integrations
				t++; //advance integration counter
			}
			t--; //prevents skipping
			snapshotLSTindices.push_back(integrationIndices); //accumulate snapshots
		}
	}
	nSnapshots = snapshotLSTindices.size();
	cout << "These integrations have been grouped into " << nSnapshots << " snapshots of at most " << snapshotTime << " seconds." << endl;
	return snapshotLSTindices;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// REPHASE AND RENORMALIZE
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//This function calculates N^-1 y, the first step in the mapmaking algorithm
void noiseInverseVarianceWeight(vector< vector<complex> >& allVisibilities, vector< vector<double> >& noiseVarianceOnEachVisibiltiy, vector<bool>& integrationsToUse){
	cout << "Now multiplying all visibilities by the inverse noise variance...[NOTE: FOR NOW, EACH VISIBILITY HAS A NOISE VARIANCE OF 1]" << endl;
	for (int t = 0; t < nLSTs; t++){
		if (integrationsToUse[t]){
			for (int b = 0; b < nBaselines; b++) allVisibilities[b][t] = allVisibilities[b][t]*(1.0/noiseVarianceOnEachVisibiltiy[b][t]);
		}
	}
}

//This function multiplies each visibility by e^(i*b.k_0) where b is the baseline vector and k_0 points to the facet center
void rephaseVisibilities(vector< vector<complex> >& allVisibilities, vector<bool>& integrationsToUse, vector<cartVec>& baselines, vector<horizPoint>& facetCenterPointings){
	cout << "Now rephasing all visibilities to the facet center..." << endl;
	cout << "[NOTE: DESPITE WHAT I ORIGINALLY THOUGHT, THE SIGN ON THE TRANSLATION SEEMS TO BE POSITIVE NOT NEGATIVE. THIS WARRENTS FURTHER INVESTIGATION.]" << endl;
	for (int t = 0; t < nLSTs; t++){
		if (integrationsToUse[t]){
			for (int b = 0; b < nBaselines; b++){
				double b_dot_k = baselines[b].dot(facetCenterPointings[t].toVec())*(k); //negative sign converts to incoming wave...TODO: or does it? It seems like it should be positive all along?!?!?
				complex complexFactor = complex(cos(b_dot_k), sin(b_dot_k)); //In A we multiply by e^-b.k so in A^t we multiply by e^b.k
				allVisibilities[b][t] = allVisibilities[b][t]*(complexFactor);
			}
		}
	}
}

void convertToTemperatureUnits(vector< vector<complex> >& allVisibilities, vector<bool>& integrationsToUse){
	cout << "Now converting the visibilities to temperature units..." << endl;
	for (int t = 0; t < nLSTs; t++){
		if (integrationsToUse[t]){
			for (int b = 0; b < nBaselines; b++) allVisibilities[b][t] = allVisibilities[b][t] * temperatureConversionFactor;
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GRIDDING VISIBILITIES TOGETHER AND FFT
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//This function figures out the unit vectors for the direction of the pixelization either in the plane of the array for the FFT or in the facet for the purposes of the primary beam calculation
//They are stored as principalAxesOfTheFacet[+RA, +Dec, +RA cross +Dec]
vector<cartVec> findFacetPrincipalAxes(double LST, vector<cartVec>& arrayPrincipalAxes, bool projectOntoArrayPlane){
	cartVec center = (equaPoint(facetRAinRad, facetDecInRad)).toHoriz(LST).toVec();
	cartVec plusDec = (equaPoint(facetRAinRad, facetDecInRad + 1e-7)).toHoriz(LST).toVec();
	cartVec plusRA = (equaPoint(facetRAinRad + 1e-7, facetDecInRad)).toHoriz(LST).toVec();

	cartVec	plusRADirection = plusRA - center;
	cartVec plusDecDirection = plusDec - center;

	if (projectOntoArrayPlane){
		plusRADirection = (plusRADirection - arrayPrincipalAxes[2]*(plusRADirection.dot(arrayPrincipalAxes[2])));
		plusDecDirection = (plusDecDirection - arrayPrincipalAxes[2]*(plusDecDirection.dot(arrayPrincipalAxes[2])));
	}
	cartVec RAcrossDecDirection = plusRADirection.cross(plusDecDirection);

	vector<cartVec> principalAxesOfTheFacet;
	principalAxesOfTheFacet.push_back(plusRADirection.normalize());
	principalAxesOfTheFacet.push_back(plusDecDirection.normalize());
	principalAxesOfTheFacet.push_back(RAcrossDecDirection.normalize());
	return principalAxesOfTheFacet;
}

//This function figures out the size of the grid we're going to need to cover the whole celestial sphere to prevent aliasing. 
//Consequently, it figures out the proper pixelization of the baseline vectors b_alpha and b_delta
//The way it accomplishes this is by going around the unit circle in degree increments and figuring out how far we'd need to go in the alpha and delta directions
vector<int> computeFullFieldSizes(int& nAlphaMax, int& nDeltaMax, double LST, vector<cartVec>& facetAxesInArrayPlane, vector<cartVec>& arrayPrincipalAxes){
	int nAlpha = 2;
	int nDelta = 2;
	cartVec toFacetCenter = ((equaPoint(facetRAinRad, facetDecInRad)).toHoriz(LST)).toVec();
	cartVec	toFacetCenterProjectionOntoArrayPlane = toFacetCenter - arrayPrincipalAxes[2]*(toFacetCenter.dot(arrayPrincipalAxes[2]));
	double alphaDotDelta = facetAxesInArrayPlane[0].dot(facetAxesInArrayPlane[1]);

	for (double theta = 0; theta < 2*pi; theta+=(pi/180.0)){
		cartVec fromProjectedCenterToEdge = arrayPrincipalAxes[0]*cos(theta) + arrayPrincipalAxes[1]*sin(theta) - toFacetCenterProjectionOntoArrayPlane;
		double dotIntoAlpha = fromProjectedCenterToEdge.dot(facetAxesInArrayPlane[0]);
		double dotIntoDelta = fromProjectedCenterToEdge.dot(facetAxesInArrayPlane[1]);
		double alphaComp = (dotIntoAlpha - alphaDotDelta * dotIntoDelta) / (1.0 - pow(alphaDotDelta,2));
		double deltaComp = (dotIntoDelta - alphaDotDelta * dotIntoAlpha) / (1.0 - pow(alphaDotDelta,2));
		if (nAlpha < ceil(fabs(alphaComp) / angResInRad)*2) nAlpha = int(ceil(fabs(alphaComp) / angResInRad)*2);
		if (nDelta < ceil(fabs(deltaComp) / angResInRad)*2) nDelta = int(ceil(fabs(deltaComp) / angResInRad)*2);
	}

	if (nAlphaMax < nAlpha) nAlphaMax = nAlpha;
	if (nDeltaMax < nDelta) nDeltaMax = nDelta;
	vector<int> fullFieldSize;
	fullFieldSize.push_back(nAlpha);
	fullFieldSize.push_back(nDelta);
	return fullFieldSize;
}

//This function converts baselines to array indices accoring to the formula j = n * Delta Theta * k * b_alpha/ (2*pi) and then grids up the visibilities and their complex conjugates
void addVisibilitiesToFacetGrid(vector< vector<complex> >& facetGrid, int nAlpha, int nDelta, vector<cartVec>& projectedBaselines, vector<cartVec>& arrayPrincipalAxes, int LSTindex, double LST, vector< vector<complex> >& allVisibilities){
	vector<cartVec> facetAxesInArrayPlane = findFacetPrincipalAxes(LST, arrayPrincipalAxes, true);
	for (int b = 0; b < nBaselines; b++){
		double b_alpha = (projectedBaselines[b].x * facetAxesInArrayPlane[0].x) + (projectedBaselines[b].y * facetAxesInArrayPlane[0].y);
		double b_delta = (projectedBaselines[b].x * facetAxesInArrayPlane[1].x) + (projectedBaselines[b].y * facetAxesInArrayPlane[1].y);
		int alphaIndex = int(round(k * angResInRad * b_alpha * nAlpha / (2*pi)));
		int deltaIndex = int(round(k * angResInRad * b_delta * nDelta / (2*pi)));
		if ((abs(alphaIndex) < nAlpha/2) && (abs(deltaIndex) < nDelta/2)){
			facetGrid[nAlpha/2 + alphaIndex][nDelta/2 + deltaIndex] = facetGrid[nAlpha/2 + alphaIndex][nDelta/2 + deltaIndex] + allVisibilities[b][LSTindex];
			facetGrid[nAlpha/2 - alphaIndex][nDelta/2 - deltaIndex] = facetGrid[nAlpha/2 - alphaIndex][nDelta/2 - deltaIndex] + allVisibilities[b][LSTindex].conj();
		}
	}
}

//This function takes the accumulated, gridded visibilities for a given snapshot and performs an FFT
vector< vector<double> > performUnormalizedShiftedFFT(vector< vector<complex> >& facetGrid, int nAlpha, int nDelta, bool forward){
	fftw_complex *in, *out;
	in = (fftw_complex*) fftw_malloc(nAlpha*nDelta * sizeof(fftw_complex));
	out = (fftw_complex*) fftw_malloc(nAlpha*nDelta * sizeof(fftw_complex));
	//Shift and populate FFTW object
	for (int a = 0; a < nAlpha; a++){
		for (int d = 0; d < nDelta; d++){
			int shiftedAlpha = (a + nAlpha/2)%nAlpha;
			int shiftedDelta = (d + nDelta/2)%nDelta;
			in[shiftedAlpha*nDelta + shiftedDelta][0] = facetGrid[a][d].re;
			in[shiftedAlpha*nDelta + shiftedDelta][1] = facetGrid[a][d].im;
		}
	}

	//Perform FFT
	fftw_plan FT;
	if (forward){ 
		FT = fftw_plan_dft_2d(nAlpha, nDelta, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	} else {
		FT = fftw_plan_dft_2d(nAlpha, nDelta, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	}
	fftw_execute(FT);

	//Shift and put back into facetMap
	vector< vector<double> > facetMap(nAlpha, vector<double>(nDelta, 0.0));
	for (int a = 0; a < nAlpha; a++){
		for (int d = 0; d < nDelta; d++){
			int shiftedAlpha = (a + nAlpha/2)%nAlpha;
			int shiftedDelta = (d + nDelta/2)%nDelta;
			facetMap[a][d] = out[shiftedAlpha*nDelta + shiftedDelta][0];
		}
	}

	fftw_destroy_plan(FT);
	fftw_free(in); fftw_free(out);
	return facetMap;
}


//This function retruns the value of the gain of the PB as a funciton of altitude and azimuth
//Unfortunately, the primary beam azimuth is stored with the XX polarization as azimuth zero and continues CCW. 
double primaryBeam(horizPoint& pointing, vector< vector<double> >& PB){
	if (gaussianPrimaryBeam){
		if (pointing.alt > 0){
			double sigma = primaryBeamFWHM/360.0*2*pi/2.355;
			return exp(-pow(pi/2 -pointing.alt,2)/2/pow(sigma,2))/sigma/sqrt(2*pi);
		} else {
			return 0.0;
		}
	}


	if (pointing.alt <= 0) return 0.0;

	double altPixel = pointing.alt * 180.0 / pi / deltaBeamAlt;
	double azPixel = fmod(-pointing.az * 180.0 / pi + xpolOrientationDegreesEastofNorth + 360.0,360.0) / deltaBeamAz;
	int altIndex1 = int(floor(altPixel));
	int altIndex2 = int(ceil(altPixel));
	int azIndex1 = int(floor(azPixel));
	int azIndex2 = int(ceil(azPixel));
	
	//Handle integer pixel values to avoid getting 0/0
	if ((altIndex1 == altIndex2) && (azIndex1 == azIndex2)) return (PB[altIndex1][azIndex1]);
	if (altIndex1 == altIndex2) return ((PB[altIndex1][azIndex2] - PB[altIndex1][azIndex1]) * (azPixel - azIndex1) + PB[altIndex1][azIndex1]);
	if (azIndex1 == azIndex2) return ((PB[altIndex2][azIndex1] - PB[altIndex1][azIndex1]) * (altPixel - altIndex1) + PB[altIndex1][azIndex1]);

	double PBresult = (PB[altIndex1][azIndex1] * (altIndex2-altPixel) * (azIndex2-azPixel));
	PBresult += (PB[altIndex2][azIndex1] * (altPixel-altIndex1) * (azIndex2-azPixel));
	PBresult += (PB[altIndex1][azIndex2] * (altIndex2-altPixel) * (azPixel-azIndex1));
	PBresult += (PB[altIndex2][azIndex2] * (altPixel-altIndex1) * (azPixel-azIndex1));
	return PBresult;
}

//This function figures out the angle in horizontal coordinates that points to each grid point center on the facet in 3D space.
//It then multiplies the image in snapshotMap by the interpolated value of the primary beam at that horiztonal pointing.
//This function also multiplies by the factor Delta A_k / (k * k_z). k_z is calculated by looking at k_x and k_y (as defined by the plane of the array).
//The value of k_z is also determined by looking at the horizontal pointing of pixel and seeing where that interesects the unit sphere.
//There is another way to approach this problem, which is to project straight down onto the plane of the array. I'm not sure that that's worse, so figuring that out is something TODO.
void multiplySnapshotByPrimaryBeamAndGeometricFactors(vector< vector<double> >& snapshotMap, vector< vector<double> >& PB, vector<cartVec>& arrayPrincipalAxes, int nAlpha, int nDelta, double LST){
	vector<cartVec> empty;
	vector<cartVec> principalAxesOfTheFacet = findFacetPrincipalAxes(LST, empty, false);
	vector<cartVec> facetAxesInArrayPlane = findFacetPrincipalAxes(LST, arrayPrincipalAxes, true);
	double pixelAreaElement = pow(k * angResInRad, 2) * fabs(facetAxesInArrayPlane[0].x * facetAxesInArrayPlane[1].y - facetAxesInArrayPlane[0].y * facetAxesInArrayPlane[1].x);
	cartVec facetCenterVector = (equaPoint(facetRAinRad, facetDecInRad)).toHoriz(LST).toVec();
	for (int a = 0; a < nAlpha; a++){
		for (int d = 0; d < nDelta; d++){
			cartVec	positionInFacet = facetCenterVector + principalAxesOfTheFacet[0]*((a - nAlpha/2)*angResInRad);
			positionInFacet = positionInFacet + principalAxesOfTheFacet[1]*((d - nDelta/2)*angResInRad);
			horizPoint gridPointDirection = positionInFacet.toHoriz();
			double k_z = k * sin(gridPointDirection.alt);
			snapshotMap[a][d] *= (primaryBeam(gridPointDirection, PB) * pixelAreaElement / (k * k_z));
		}
	}
}

// This function addes a snapshotMap, which might be smaller thant the largest possible map, to the coaddedMap, which combines all properly weighted snapshots.
void addSnapshotMap(vector< vector<double> >& snapshotMap, vector< vector<double> >& coaddedMap, int nAlpha, int nAlphaMax, int nDelta, int nDeltaMax){
	for (int a = 0; a < nAlpha; a++){
		for (int d = 0; d < nDelta; d++){
			coaddedMap[a - nAlpha/2 + nAlphaMax/2][d - nDelta/2 + nDeltaMax/2] += snapshotMap[a][d];
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// POINT SPREAD FUNCTION CALCULATION
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//This function calculates the inverse noise variance for the gridded visibilities.
//Since variances add for the sum of uncorrelated random variables, I'm simply computing the inverse variance by adding together the variances on each visibility by assigning them to cells.
vector< vector<double> > calculateNinv(int nAlpha, int nDelta, vector<cartVec>& projectedBaselines, vector< vector<double> >& noiseVarianceOnEachVisibiltiy, vector<cartVec>& arrayPrincipalAxes, vector<int>& LSTindices, vector<double>& LSTs){
	vector< vector<double> > Ninv(nAlpha, vector<double>(nDelta, 0.0));
	for (int i = 0; i < LSTindices.size(); i++){
		vector<cartVec> facetAxesInArrayPlane = findFacetPrincipalAxes(LSTs[LSTindices[i]], arrayPrincipalAxes, true);
		for (int b = 0; b < nBaselines; b++){
			double b_alpha = (projectedBaselines[b].x * facetAxesInArrayPlane[0].x) + (projectedBaselines[b].y * facetAxesInArrayPlane[0].y);
			double b_delta = (projectedBaselines[b].x * facetAxesInArrayPlane[1].x) + (projectedBaselines[b].y * facetAxesInArrayPlane[1].y);
			int alphaIndex = int(round(k * angResInRad * b_alpha * nAlpha / (2*pi)));
			int deltaIndex = int(round(k * angResInRad * b_delta * nDelta / (2*pi)));
			if ((abs(alphaIndex) < nAlpha/2) && (abs(deltaIndex) < nDelta/2)){
				Ninv[nAlpha/2 + alphaIndex][nDelta/2 + deltaIndex] = Ninv[nAlpha/2 + alphaIndex][nDelta/2 + deltaIndex] + 1.0/noiseVarianceOnEachVisibiltiy[b][LSTindices[i]];
				Ninv[nAlpha/2 - alphaIndex][nDelta/2 - deltaIndex] = Ninv[nAlpha/2 - alphaIndex][nDelta/2 - deltaIndex] + 1.0/noiseVarianceOnEachVisibiltiy[b][LSTindices[i]];
			}
		}
	}

	return Ninv;
}

//This function looks at Ninv to figure out which values of a and d correspond to baselines that have been sampled. Then the baseline-facing part of A^t can be treated as a sparse matrix, saving memory
vector<coord> listAllSampledBaselines(int nAlpha, int nDelta, vector< vector<double> >& Ninv){
	vector<coord> sampledBaselines;
	for (int a = 0; a < nAlpha; a++){
		for (int d = 0; d < nDelta; d++){
			if (Ninv[a][d] > 0){
				coord sampled(a,d);
				sampledBaselines.push_back(sampled);
			}
		}
	}
	return sampledBaselines;
}

//This function computes A^t matrix for each snapshot. The action of this matrix is to convert gridded visibilities into a dirty map.
//The matrix wraps faster over Dec faster than RA and over y faster than x in the array
vector< vector<complex> > calculateKAtranspose(double LST, int nAlpha, int nDelta, vector<cartVec>& arrayPrincipalAxes, vector<cartVec>& projectedBaselines, vector< vector<double> >& PB, horizPoint facetCenter, vector<coord>& sampledBaselines){

	//Calculate the real space diagonal part of Atranspose
	vector< vector<double> > realSpaceDiagonalPart(nAlphaCropped, vector<double>(nDeltaCropped, temperatureConversionFactor));
	multiplySnapshotByPrimaryBeamAndGeometricFactors(realSpaceDiagonalPart, PB, arrayPrincipalAxes, nAlphaCropped, nDeltaCropped, LST); //takes advantage of previous code
	
	//Calculate a 2D DFT matrix
	int nSampled = sampledBaselines.size();
	vector< vector<complex> > Atranspose(nAlphaCropped*nDeltaCropped, vector<complex>(nSampled, complex(0,0)));
	for (int a = 0; a < nAlphaCropped; a++){
		for (int d = 0; d < nDeltaCropped; d++){
			for (int n = 0; n < nSampled; n++){
				int a2 = sampledBaselines[n].a;
				int d2 = sampledBaselines[n].d;
				double argument = (2*pi/nAlpha)*(a - nAlphaCropped/2)*(a2 - nAlpha/2) + (2*pi/nDelta)*(d - nDeltaCropped/2)*(d2 - nDelta/2);
				Atranspose[d + nDeltaCropped*a][n].re = cos(argument);
				Atranspose[d + nDeltaCropped*a][n].im = sin(argument);
			}
		}
	}
	//Calculate the baseline space diaognal part of Atranspose
	vector<cartVec> facetAxesInArrayPlane = findFacetPrincipalAxes(LST, arrayPrincipalAxes, true);
	vector<complex> diagonalPartBaselineSpace(nSampled, complex(0,0));
	cartVec	vecToFacetCenter = facetCenter.toVec() * k;
	for (int n = 0; n < nSampled; n++){
		int a = sampledBaselines[n].a;
		int d = sampledBaselines[n].d;
		double b_alpha = (nAlpha/2 - a) * 2 * pi / (nAlpha * angResInRad * k);
		double b_delta = (nDelta/2 - d) * 2 * pi / (nDelta * angResInRad * k);
		double divisor = (facetAxesInArrayPlane[0].y * facetAxesInArrayPlane[1].x) - (facetAxesInArrayPlane[0].x * facetAxesInArrayPlane[1].y);
		double b_x = (b_delta * facetAxesInArrayPlane[0].y - b_alpha * facetAxesInArrayPlane[1].y) / divisor;
		double b_y = (b_alpha * facetAxesInArrayPlane[1].x - b_delta * facetAxesInArrayPlane[0].x) / divisor;
		cartVec pixelizedBaseline = (arrayPrincipalAxes[0]*b_x) + (arrayPrincipalAxes[1]*b_y);
		double b_dot_k = pixelizedBaseline.dot(vecToFacetCenter);
		diagonalPartBaselineSpace[n].re = cos(-b_dot_k); //TODO: ALTHOUGH THIS APPEARS TO BE THE RIGHT SIGN, I'M NOT SURE I CAN JUSTIFY IT.
		diagonalPartBaselineSpace[n].im = sin(-b_dot_k);
	}

	//Multiply the three together
	for (int a = 0; a < nAlphaCropped; a++){
		for (int d = 0; d < nDeltaCropped; d++){
			for (int n = 0; n < nSampled; n++){
				int a2 = sampledBaselines[n].a;
				int d2 = sampledBaselines[n].d;
				Atranspose[d + nDeltaCropped*a][n] = Atranspose[d + nDeltaCropped*a][n] * diagonalPartBaselineSpace[n] * realSpaceDiagonalPart[a][d];
			}
		}
	}

	return Atranspose;
}



//This function computes KAtNinvAKt for each snapshot and then adds them to the existing sum to get an overall PSF
void addSnapshotPSF(vector< vector<double> >& PSF, vector< vector<complex> >& KAtranspose, vector< vector<double> >& Ninv, int nAlpha, int nDelta, vector<coord>& sampledBaselines){
	vector< vector<double> > snapshotPSF(nAlphaFacet*nDeltaFacet, vector<double>(nAlphaCropped*nDeltaCropped, 0.0));
	for (int n = 0; n < sampledBaselines.size(); n++){
		int a = sampledBaselines[n].a;
		int d = sampledBaselines[n].d;
		for (int af = 0; af < nAlphaFacet; af++){
			for (int df = 0; df < nDeltaFacet; df++){
				int facetIndexInExtendedMap = (PSFextensionBeyondFacetFactor-1)/2 * nDeltaFacet + df + nDeltaCropped * ((PSFextensionBeyondFacetFactor-1)/2 * nAlphaFacet + af);
				for (int j = 0; j < nAlphaCropped*nDeltaCropped; j++){
					PSF[df+nDeltaFacet*af][j] = PSF[df+nDeltaFacet*af][j] + (KAtranspose[facetIndexInExtendedMap][n] * KAtranspose[j][n].conj() * Ninv[a][d]).re; //TODO: check imaginary component 
				}
			}
		}
	}
}

//This function normalizes the PSF to 1 at the desired point probed and it computes the diagonal normalization matrix needed to do that.
vector<double> computeNormalizationAndNormalizePSF(vector< vector<double> >& PSF){
	cout << "Now normalizing the PSF..." << endl;
	vector<double> normalizationMatrix(nAlphaFacet*nDeltaFacet,0.0);
	vector< vector<double> > PSFtemp = PSF;
	for (int af = 0; af < nAlphaFacet; af++){
		for (int df = 0; df < nDeltaFacet; df++){
			if (PSFpeaksAtOne) {
				double maxPSFValueInThisRow = 0.0;
				for (int j = 0; j < nAlphaCropped*nDeltaCropped; j++){
					if (PSF[df+nDeltaFacet*af][j] > maxPSFValueInThisRow) maxPSFValueInThisRow = PSF[df+nDeltaFacet*af][j];
				}
				normalizationMatrix[df+nDeltaFacet*af] = 1.0/maxPSFValueInThisRow;
			} else {
				int facetIndexInExtendedMap = (PSFextensionBeyondFacetFactor-1)/2 * nDeltaFacet + df + nDeltaCropped * ((PSFextensionBeyondFacetFactor-1)/2 * nAlphaFacet + af);
				normalizationMatrix[df+nDeltaFacet*af] = 1.0/PSFtemp[df+nDeltaFacet*af][facetIndexInExtendedMap];	
			}
			for (int j = 0; j < nAlphaCropped*nDeltaCropped; j++) PSF[df+nDeltaFacet*af][j] = normalizationMatrix[df+nDeltaFacet*af]*PSFtemp[df+nDeltaFacet*af][j];
		}
	}
	return normalizationMatrix;
}

//This function renormalizes the dirty map so that it is the true map convolved with the PSF (on average).
vector< vector<double> > renormalizeMap(vector<double>& normalizationMatrix, vector< vector<double> >& coaddedMap, int nAlphaMax, int nDeltaMax){
	vector< vector<double> > renormalizedMap(nAlphaCropped, vector<double>(nDeltaCropped,0.0));
	for (int a = 0; a < nAlphaFacet; a++){
		for (int d = 0; d < nDeltaFacet; d++){
			renormalizedMap[a][d] = coaddedMap[nAlphaMax/2 + (a-nAlphaFacet/2)][nDeltaMax/2 + (d-nDeltaFacet/2)] * normalizationMatrix[d + nDeltaFacet*a];
		}
	}
	return renormalizedMap;
}

vector< vector<double> > computeNoiseCovariance(vector< vector<double> >& PSF, vector<double>& normalizationMatrix){
	vector< vector<double> > noiseCovariance(nAlphaFacet*nDeltaFacet, vector<double>(nAlphaFacet*nDeltaFacet,0.0));
	for (int i = 0; i < nAlphaFacet*nDeltaFacet; i++){
		for (int af = 0; af < nAlphaFacet; af++){
			for (int df = 0; df < nDeltaFacet; df++){
				int facetIndexInExtendedMap = (PSFextensionBeyondFacetFactor-1)/2 * nDeltaFacet + df + nDeltaCropped * ((PSFextensionBeyondFacetFactor-1)/2 * nAlphaFacet + af);
				noiseCovariance[i][df + nDeltaFacet*af] = PSF[i][facetIndexInExtendedMap]*normalizationMatrix[df + nDeltaFacet*af];
			}
		}
	}
	return noiseCovariance;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OUTPUT FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void exportMatrix(vector< vector<double> >& matrix, int dim1, int dim2, string filename){
	cout << "Now saving " << filename << " as a " << dim1 << " by " << dim2 << " matrix..." << endl;
	ofstream outfile;
	outfile.precision(16);	
	outfile.open(filename.c_str(), ios::trunc);		
	for (int i = 0; i < dim1; i++){
		for (int j = 0; j < dim2; j++){
			outfile << matrix[i][j] << " ";
		}
		outfile << endl;
	} 
	outfile.close();
}

void exportVector(vector<double>& vec, int dim, string filename){
	cout << "Now saving " << filename << " as a size " << dim << " vector..." << endl;
	ofstream outfile;
	outfile.precision(16);	
	outfile.open(filename.c_str(), ios::trunc);		
	for (int i = 0; i < dim; i++) outfile << vec[i] << endl;
	outfile.close();
}


void exportComplexMatrix(vector< vector<complex> >& matrix, int dim1, int dim2, string filenameReal, string filenameImag){
	cout << "Now saving " << filenameReal << " as the " << dim1 << " by " << dim2 << " as the real matrix..." << endl;
	ofstream outfile;
	outfile.precision(16);	
	outfile.open(filenameReal.c_str(), ios::trunc);		
	for (int j = 0; j < dim2; j++){
		for (int i = 0; i < dim1; i++){
			outfile << matrix[i][j].re << " ";
		}
		outfile << endl;
	} 
	outfile.close();

	cout << "Now saving " << filenameImag << " as the " << dim1 << " by " << dim2 << " as the imaginary matrix..." << endl;
	ofstream outfile2;
	outfile2.precision(16);	
	outfile2.open(filenameImag.c_str(), ios::trunc);		
	for (int j = 0; j < dim2; j++){
		for (int i = 0; i < dim1; i++){
			outfile2 << matrix[i][j].im << " ";
		}
		outfile2 << endl;
	} 
	outfile2.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN FUNCTION
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]){
	
	cout << endl << "Running Faceted Mapmaking..." << endl << endl;
	//Load relevant quanties and data
	if (argc == 2) specsFilename = argv[1];
	loadSpecs();
	vector<int> baselineRedundancy;
	vector<cartVec> baselines = loadBaselines(baselineRedundancy);
	vector<cartVec> arrayPrincipalAxes = loadArrayPrincipalAxes();
	vector<cartVec> projectedBaselines = calculateProjectedBaselines(baselines,arrayPrincipalAxes);
	vector<double> LSTs; //LSTs of each integration
	vector< vector<complex> > allVisibilities = loadVisibilities(baselines, LSTs);
	vector< vector<double> > noiseVarianceOnEachVisibiltiy = loadNoiseVarianceOnEachVisibiltiy(baselines, LSTs, baselineRedundancy);
	vector<horizPoint> PBpointings = loadPBPointings(LSTs);
	vector< vector<double> > discretizedPrimaryBeam = loadDiscretizedPrimaryBeam();
	
	//Geometric calculations 
	vector<horizPoint> facetCenterPointings = computeFacetCenterPointings(LSTs);
	vector<bool> integrationsToUse = determineIntegrationsToUse(PBpointings, facetCenterPointings);
	if (nIntegrationsToUse == 0) return 1;
	vector< vector<int> > snapshotLSTindices = assignLSTsToSnapshots(integrationsToUse);

	//Rephase and renormalize
	noiseInverseVarianceWeight(allVisibilities, noiseVarianceOnEachVisibiltiy, integrationsToUse);
	rephaseVisibilities(allVisibilities, integrationsToUse, baselines, facetCenterPointings);
	convertToTemperatureUnits(allVisibilities, integrationsToUse);
	
	//Loop over all snapshots to figure out the pixelization
	vector< vector<int> > fullFieldSizes;
	int nAlphaMax = 2;
	int nDeltaMax = 2;
	cout << "Now figuring out the full size of the facet for all snapshots..." << endl;
	for (int n = 0; n < nSnapshots; n++){
		double snapshotCentralLST = LSTs[snapshotLSTindices[n][int(round(snapshotLSTindices[n].size()/2.0-.5))]];
		vector<cartVec> principalAxesOfTheFacetInThePlaneOfTheArray = findFacetPrincipalAxes(snapshotCentralLST, arrayPrincipalAxes, true);
		fullFieldSizes.push_back(computeFullFieldSizes(nAlphaMax, nDeltaMax, snapshotCentralLST, principalAxesOfTheFacetInThePlaneOfTheArray, arrayPrincipalAxes));
	}

	//Grid and FFT each snapshot
	vector< vector<double> > coaddedMap(nAlphaMax, vector<double>(nDeltaMax, 0.0));
	cout << "Now calculating all snapshot maps..." << endl;
	for (int n = 0; n < nSnapshots; n++){		
		cout << " " << floor(100.0 * n / nSnapshots) << "% done. \r" << std::flush;
		int nAlpha = fullFieldSizes[n][0];
		int nDelta = fullFieldSizes[n][1];
		double snapshotCentralLST = LSTs[snapshotLSTindices[n][int(round(snapshotLSTindices[n].size()/2.0-.5))]];
		vector< vector<complex> > facetGrid(nAlpha, vector<complex>(nDelta, complex(0,0)));
		for (int i = 0; i < snapshotLSTindices[n].size(); i++){
			addVisibilitiesToFacetGrid(facetGrid, nAlpha, nDelta, projectedBaselines, arrayPrincipalAxes, snapshotLSTindices[n][i], LSTs[snapshotLSTindices[n][i]], allVisibilities);
		}
		vector< vector<double> > snapshotMap = performUnormalizedShiftedFFT(facetGrid, nAlpha, nDelta, false);
		multiplySnapshotByPrimaryBeamAndGeometricFactors(snapshotMap, discretizedPrimaryBeam, arrayPrincipalAxes, nAlpha, nDelta, snapshotCentralLST);
		addSnapshotMap(snapshotMap, coaddedMap, nAlpha, nAlphaMax, nDelta, nDeltaMax);
	}
	cout << "Done.                  " << endl;  
	exportMatrix(coaddedMap, nAlphaMax, nDeltaMax, unnormalizedFullMapFilename);

	//Calculate the cropped PSF
	cout << "Now calculating the position-dependent point spread function..." << endl;
	vector< vector<double> > PSF(nAlphaFacet*nDeltaFacet, vector<double>(nAlphaCropped*nDeltaCropped,0));
	for (int n = 0; n < nSnapshots; n++){		
		int nAlpha = fullFieldSizes[n][0];
		int nDelta = fullFieldSizes[n][1];
		cout << " " << floor(100.0 * n / nSnapshots) << "% done. \r" << std::flush;
		int snapshotCentralLSTindex = snapshotLSTindices[n][int(round(snapshotLSTindices[n].size()/2.0-.5))];
		vector< vector<double> > Ninv = calculateNinv(nAlpha, nDelta, projectedBaselines, noiseVarianceOnEachVisibiltiy, arrayPrincipalAxes, snapshotLSTindices[n], LSTs);
		vector<coord> sampledBaselines = listAllSampledBaselines(nAlpha, nDelta, Ninv);
		vector< vector<complex> > KAtranspose = calculateKAtranspose(LSTs[snapshotCentralLSTindex], nAlpha, nDelta, arrayPrincipalAxes, projectedBaselines, discretizedPrimaryBeam, facetCenterPointings[snapshotCentralLSTindex], sampledBaselines);
		addSnapshotPSF(PSF, KAtranspose, Ninv, nAlpha, nDelta, sampledBaselines); 
	}
	cout << "Done.                  " << endl;  
	vector<double> normalizationMatrix = computeNormalizationAndNormalizePSF(PSF);
	exportVector(normalizationMatrix,nAlphaFacet*nDeltaFacet,DmatrixFilename);

	//Compute and save the final data products
	exportMatrix(PSF, nAlphaFacet*nDeltaFacet, nAlphaCropped*nDeltaCropped, PSFfilename);
	vector< vector<double> > renormalizedMap = renormalizeMap(normalizationMatrix, coaddedMap, nAlphaMax, nDeltaMax);
	exportMatrix(renormalizedMap, nAlphaFacet, nDeltaFacet, finalMapFilename);
	vector< vector<double> > noiseCovariance = computeNoiseCovariance(PSF, normalizationMatrix);
	exportMatrix(noiseCovariance, nAlphaFacet*nDeltaFacet, nAlphaFacet*nDeltaFacet, noiseCovFilename);


	cout << "Done making a map. " << endl << endl;
	return 0;
}
