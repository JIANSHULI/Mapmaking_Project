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
double south, east, up, freq, startingLST, LST, duration, timeStep, latitude, longitude, xpolOrientationDegreesEastofNorth, UTCtime;
string specsFile = "Specifications.txt";
string componentsFile = "components.dat";
int componentsToFit = 9;
string FITS_directory, alm_directory, GSMformat1, GSMformat2, polarization, dateString, timeString, resultsFolder;
int startingPower, endingPower;
int principalComps = 3; //Hardcoded...won't work otherwise
bool beamIsOne = false;

//Loads the relevant information about the GSM files and the observation into global variables
void loadSpecs(){
	fstream infile(specsFile.c_str(),fstream::in);
	string dummy;
	infile >> dummy >> FITS_directory;
	infile >> dummy >> GSMformat1;
	infile >> dummy >> GSMformat2;
	infile >> dummy >> startingPower;
	infile >> dummy >> endingPower;
	infile >> dummy >> latitude;
	infile >> dummy >> longitude;
	infile >> dummy >> xpolOrientationDegreesEastofNorth;
	infile >> dummy >> resultsFolder;
	infile.close();
}

//Determine which map satisfies the relationship k*d > nSides; nSide = 8 disallowed
int pickResolution(double wavenumber){
	for (int res = startingPower*2; res <= endingPower; res*=2){
		if (wavenumber*sqrt(pow(south,2)+pow(east,2)) < res){
			cout << "GSM" << res << " selected." << endl;
			return res;
		}
	}
	cout << endl << "WARNING: The maximum resolution GSM available is not high enough resolution." << endl << endl;
	return endingPower;
}

//Loads in all 3 principal components of the GSM and return them as a vector.
vector< Healpix_Map<double> > loadGSMpcomps(int res){
	vector< Healpix_Map<double> > GSMpcomps;
	for (int p = 1; p <= principalComps; p++){
		stringstream ss;
		ss << "./" << FITS_directory << GSMformat1 << p << GSMformat2 << res;
		string GSMfile = ss.str();
		cout << "Now loading " << GSMfile << endl;
		Healpix_Map<double> map; 
		read_Healpix_map_from_fits(GSMfile,map); 
		GSMpcomps.push_back(map);
	}
	return GSMpcomps;
}

//Load the components data and properly reweight and rescale the GSM principal componenets
Healpix_Map<double> computeGSM(vector< Healpix_Map<double> >& GSMpcomps, int res){
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

/* START SHANA'S BEAM CODE */
	//Jeff's helper code from Shana's loading code
	string exec(string input_command) {//return result of stdout of a system command such as "python BLAH_BLAH". Won't catch stderr. Direct copy from internet http://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c
		//cout << " test-2 " + input_command;
		char* cmd = (char*)input_command.c_str();
		//cout << " test-1 " + input_command;
		//cout << cmd << endl;
		FILE* pipe = popen(cmd, "r");
		if (!pipe) return "ERROR";
		char buffer[1024];
		std::string result = "";
		//cout << " test0 " + input_command;
		while(!feof(pipe)) {
			if(fgets(buffer, 1024, pipe) != NULL){
				//cout << buffer << " ";			
				result += buffer;
			}	
		}
		pclose(pipe);
		if ( result.size() > 0 )
			result.erase(result.end()-1,result.end());
		//cout << " test1 " + input_command;
		return result;
	}
	//Jeff's helper code from Shana's loading code
	vector<double> strtovf(string in){
		//cout << "DEBUG " << in << endl;
		stringstream stream(in);
		vector<double> output;
		double tmpf;
		while (stream.good()){
			stream >> tmpf;
			//cout << tmpf << endl;
			output.push_back(tmpf);
		}
		//output.pop_back(); //sometimes get a 0 in end..tricky...
		return output;
	}
	//Jeff's helper code from Shana's loading code
	int cmdCountLine(const char* inputfilename){
		string command(inputfilename);
		command = "wc " + command;
		vector<double> results = strtovf(exec(command));
		return (int)(results[0]);
	}
	//Shana's loading function for the primary beam of an antenna polarization at a given frequency
	vector<vector<vector<double> > > makeDataVectorFromFandPol(double f, char* pol) {
	  ifstream dataRead;
	  int lineNum = 0;
	  int numData = 0; //right now these are not used
	  vector<double> dataf(numData);//right now these are not used
	  int numFreqFiles = 0; 
	  vector<double> columns(8,0.0);

	  // open data for X or Y polarization
	  int result=  strcmp(pol,"x");
	  if(result<1){
	    //cout << pol << " X " << endl;
	    dataRead.open("mwa_dataX");
	    lineNum = cmdCountLine("mwa_dataX");
	  }
	  if(result>0){
	    //cout << pol << " Y " << endl;
	    dataRead.open("mwa_dataY");
	    lineNum = cmdCountLine("mwa_dataY");
	  }
	  // create data block for X or Y polarization
	  numData=lineNum*8;
	  numFreqFiles=lineNum/8100; 
	  //cout << numFreqFiles << " freqs.line# " << lineNum<<endl; 
	  vector<vector<double> > dataForOneFreq(lineNum, columns);
	  vector<vector<double> > dataForOneFreqBlock(8100,columns);
	  vector<vector<vector<double> > > dataForFile (numFreqFiles, dataForOneFreq);
	  // read data into vector
	  for(int nf = 0; nf<numFreqFiles; nf++){
	    for(int row = 0; row < 8100; row++){
	      for(int col = 0; col < 8; col++){
		dataRead >> dataForFile[nf][row][col];
		//cout << dataForFile[nf][row][col] << "\t" <<nf << "\t" << row << "\t"<< col<< endl; 
	      }
	    }
	  }
	  int dataOutN = 0;
	  int offsetFreq=0; int n = 0; int dataBlockSize=8100;
	  if (int(round(f)) == f && (int(round(f))-100) % 10 == 0){ // if freq data already known
	    dataOutN = 1;
	    offsetFreq=((f-100)/10)-1;
	    //cout << "Offset is " << offsetFreq<<endl;
	  }
	  else{ // if interpolating gain between frequencies:
	    dataOutN=2;
	    // set offset to be lower frequency block
	    // and copy the next frequency block as well
	    for(n=11;n<20;n++){
	      if((f-(n*10))<10){
	    	offsetFreq=(n-11);//*dataBlockSize;
	    	//cout << "Offset for doubleBlock is " << offsetFreq<<endl;
	    	break;
	      }
	    }
	  }
	  vector<vector<vector<double> > > dataOut(dataOutN,dataForOneFreqBlock);
	  // read data into outfile
	  if(dataOutN==1){
	    for(int row = 0; row < 8100; row++){
	      for(int col = 0; col < 8; col++){
		dataOut[0][row][col] = dataForFile[offsetFreq][row][col];
		//cout << dataOut[0][row][col] << "\t" << row << "\t"<< col<< endl; 	
	      }
	    }
	  }
	  else if(dataOutN==2){
	    for(int nf = 0; nf<2; nf++){
	      for(int row = 0; row < 8100; row++){
		for(int col = 0; col < 8; col++){
		  dataOut[nf][row][col] = dataForFile[nf+offsetFreq][row][col];
		  //cout << dataOut[nf][row][col] << "\t" <<nf << "\t" << row << "\t"<< col<< endl; 
		}
	      }
	    }
	  }
	  
	  dataRead.close();
	  return dataOut;
	}
	//Shana's beam interpolation function (somewhat confusingly, theta is azimuth, phi is altitude, both in degrees)
	double interpolateBeam(vector<vector<vector<double> > >& beamData, double f, double theta, double phi){
	  double returnGain = 0.0,newGain1 = 0.0, newGain2 = 0.0;
	  int nf = beamData.size(); //cout << nf << " is #freq sets" << endl;
	  int lineNum = beamData[0].size(); //cout << lineNum<< "=lineNum"<< endl;
	  int col = beamData[0][0].size(); //cout << col<<" is #cols"<<endl;

	  //syphon off alt az and gain columns into faster linear arrays

	  // increments for boundary value direct element lookup.
	  double deltaX = 360.0/179.0; // cout<<"deltaX = "<< deltaX<<endl;
	  double deltaY = 90.0/44.0; // cout<<"deltaY = "<< deltaY<<endl;
	  // get list of unique x and y values for sorting new points.
	  double xUnique[180], yUnique[45];
	  int xInc = 45;

	  int xNewIndexAlt = 0, yNewIndexAlt = 0; 
	  double xUpperAlt = 0.0, xLowerAlt = 0.0, yUpperAlt = 0.0, yLowerAlt = 0.0;
	  // use xnew and ynew to find boundary indices
	  xNewIndexAlt = int(theta/deltaX); //cout<<"xNewIndexAlt "<<xNewIndexAlt << endl;
	  yNewIndexAlt = int(phi/deltaY); //cout<<"yNewIndexAlt "<<yNewIndexAlt << endl;
	  xLowerAlt=beamData[0][xInc*xNewIndexAlt][0]; xUpperAlt = beamData[0][xInc*(xNewIndexAlt+1)][0];
	  yLowerAlt=beamData[0][yNewIndexAlt][1]; yUpperAlt = beamData[0][yNewIndexAlt+1][1];

	  // Apply boundary values
	  double x1 = beamData[0][(45*(xNewIndexAlt))+(yNewIndexAlt)][0];
	  double x2 = beamData[0][(45*(xNewIndexAlt+1))+(yNewIndexAlt+1)][0];
	  double y1 = beamData[0][(45*(xNewIndexAlt))+(yNewIndexAlt)][1];
	  double y2 = beamData[0][(45*(xNewIndexAlt+1))+(yNewIndexAlt+1)][1];

	  double Q11 = beamData[0][(45*(xNewIndexAlt))+(yNewIndexAlt)][6];
	  double Q21 = beamData[0][45*(xNewIndexAlt+1)+(yNewIndexAlt)][6];
	  double Q12 = beamData[0][45*(xNewIndexAlt)+(yNewIndexAlt+1)][6];
	  double Q22 = beamData[0][(45*(xNewIndexAlt+1))+(yNewIndexAlt+1)][6];

	  double x = theta, y=phi;
	  newGain1 = (1/((x2-x1)*(y2-y1)))*(Q11*(x2-x)*(y2-y)+Q21*(x-x1)*(y2-y)+Q12*(x2-x)*(y-y1)+Q22*(x-x1)*(y-y1));
	  returnGain = newGain1;
	  
	  if(nf>1){
	    x1 = beamData[1][(45*(xNewIndexAlt))+(yNewIndexAlt)][0];
	    x2 = beamData[1][(45*(xNewIndexAlt+1))+(yNewIndexAlt+1)][0];
	    y1 = beamData[1][(45*(xNewIndexAlt))+(yNewIndexAlt)][1];
	    y2 = beamData[1][(45*(xNewIndexAlt+1))+(yNewIndexAlt+1)][1];
  
	    Q11 = beamData[1][(45*(xNewIndexAlt))+(yNewIndexAlt)][6];
	    Q21 = beamData[1][45*(xNewIndexAlt+1)+(yNewIndexAlt)][6];
	    Q12 = beamData[1][45*(xNewIndexAlt)+(yNewIndexAlt+1)][6];
	    Q22 = beamData[1][(45*(xNewIndexAlt+1))+(yNewIndexAlt+1)][6];
	    newGain2 = (1/((x2-x1)*(y2-y1)))*(Q11*(x2-x)*(y2-y)+Q21*(x-x1)*(y2-y)+Q12*(x2-x)*(y-y1)+Q22*(x-x1)*(y-y1));

	  }
	  //interpolate between two gain values if needed
	  double freqLow = 0, freqHigh = 0;
	  double interpolatedGain = 0.0;
	  if(nf>1){
	    //returnGain=0;cout << "return gain value at begin "<<returnGain<<endl;
	    //cout << "interpolate between two gain values"<<endl;
	    // get bounding frequency values
	    double remainder = f - floor(f) + (int(floor(f))-100)% 10; 
	    //cout <<remainder<<"=multiple"<<endl;
	    freqLow = freq-remainder; //cout << freqLow<<endl;
	    freqHigh = freqLow+10; //cout << freqHigh<<endl;cout <<freq<<endl;
	    // get interpolated value with freqLow,freqHigh,newGain,newGain2
	    interpolatedGain = newGain1+((newGain2-newGain1)*((f-freqLow)/(freqHigh-freqLow)));
	    //cout<<"return interpolated gain is "<<interpolatedGain<<endl;
	    returnGain = interpolatedGain;
	  }
	  return returnGain;
	}
/* END SHANA'S BEAM CODE */

//Converts vector direction of x (toward Gal Center), y, z (toward NGP) to (altitude, azimuth east from north) at the specified location and local siderial time
pointing convertGalToHoriz(pointing& galPointing, bool verbose){

	double b = pi/2.0 - galPointing.theta;
	double l = galPointing.phi;
	double pole_ra = 2.0*pi/360.0*192.859508;
    double pole_dec = 2.0*pi/360.0*27.128336;
    double posangle = 2.0*pi/360.0*(122.932-90.0);

    double ra = atan2((cos(b)*cos(l-posangle)), (sin(b)*cos(pole_dec) - cos(b)*sin(pole_dec)*sin(l-posangle))) + pole_ra;
    double dec = asin(cos(b)*cos(pole_dec)*sin(l-posangle) + sin(b)*sin(pole_dec));

    if (verbose){
    	cout << "LST: " << LST << endl;
    	cout << "l: " << 360.0/2/pi*l << endl;
    	cout << "b: " << 360.0/2/pi*b << endl;
    	cout << "ra: " << 360.0/2/pi*ra << endl;
    	cout << "dec: " << 360.0/2/pi*dec << endl;
    }

    double lha = 2.0*pi/24.0*(LST - ra*24.0/2.0/pi);
    double latInRad = 2.0*pi/360.0*latitude;
    pointing horizPointing;
    
    // performing conversion to horizontal coordinate system
    //horizPointing.phi = atan2(sin(lha), cos(lha) * sin(latInRad) - tan(dec) * cos(latInRad)); //Azimuth
    horizPointing.phi = atan2(sin(lha) * cos(dec), cos(lha) * cos(dec) * sin(latInRad) - sin(dec) * cos(latInRad)) + pi; //Azimuth
    horizPointing.theta = asin(sin(latInRad) * sin(dec) + cos(latInRad) * cos(dec) * cos(lha)); //Altitude

	if (verbose){
    	cout << "az: " << 360.0/2/pi*horizPointing.phi << endl;
    	cout << "alt: " << 360.0/2/pi*horizPointing.theta << endl;
    }

	return horizPointing;
}

//Returns the value of the primary beam for a given pointing in horizontal coordiantes
double primaryBeam(pointing& horizPointing, vector<vector<vector<double> > >& beamDataX, vector<vector<vector<double> > >& beamDataY){
	//if (beamIsOne) return 1.0;
	if (horizPointing.theta > 0){
		//return pow(sin(horizPointing.theta),2);
		//Right now, 90 degrees gives the maximum beam for xPol, meaning that xPol is oriented along the 0-180 axis
		double beamTheta = -360/2/pi*horizPointing.phi + xpolOrientationDegreesEastofNorth;
		if (beamTheta < 0) beamTheta += 360.0;
		if (beamTheta > 360) beamTheta -= 360.0;
		double beamPhi = 360/2/pi*horizPointing.theta;
		if (beamPhi == 90) beamPhi += .001;
		if (polarization.compare("xx") == 0){
			return interpolateBeam(beamDataX, freq, beamTheta, beamPhi);
		} else if ((polarization.compare("xy") == 0) || polarization.compare("yx") == 0) {
			return sqrt(interpolateBeam(beamDataY, freq, beamTheta, beamPhi) * interpolateBeam(beamDataX, freq, beamTheta, beamPhi));
		} else if (polarization.compare("yy") == 0){
			return interpolateBeam(beamDataY, freq, beamTheta, beamPhi);
		} else {
			return 0.0;
		}
	} else {
		return 0.0;
	}
}

//Returns the dot product between the baseline vector and the pointing vector in horizontal coordinates
double baselineDotProduct(pointing& horizPointing){
	double pointUp = sin(horizPointing.theta);
	double pointEast = cos(horizPointing.theta) * sin(horizPointing.phi);
	double pointSouth = cos(horizPointing.theta) * (-1.0) * cos(horizPointing.phi);
	return -pointUp*up - pointEast*east - pointSouth*south;
}

//Computes complex visibilities from the GSM
xcomplex<double> computeVisibility(Healpix_Map<double>& GSM, vector<vector<vector<double> > >& beamDataX, vector<vector<vector<double> > >& beamDataY, int res){
	Healpix_Map<double> realIntegrand = GSM;
	Healpix_Map<double> imagIntegrand = GSM;
	int N = 12*res*res;
	//multiply by the beam factor and the baseline exponential by figuring out the horizontal pointing of each pixel
	//Trafo transformation(2000, 2000, Galactic, Equatorial);
	for (int n = 0; n < N; n++){
		//vec3 thisVec = GSM.pix2vec(n);
		pointing thisPoint = GSM.pix2ang(n);
		//if (GSM[n] != 0) cout << n << endl << "GalTheta: " << 360.0/2/pi*GSM.pix2ang(n).theta << endl << "GalPhi: " << 360.0/2/pi*GSM.pix2ang(n).phi << endl << endl;
		pointing horizPointing = convertGalToHoriz(thisPoint, false);
		//if (GSM[n] != 0) cout << n << endl << "theta: " << 360.0/2/pi*horizPointing.theta << endl << "phi: " << 360.0/2/pi*horizPointing.phi << endl << endl;
		//if (GSM[n] != 0) horizPointing = convertGalToHoriz(thisPoint, true);
		double beam = primaryBeam(horizPointing, beamDataX, beamDataY);
		
		double baselineDotkHat = baselineDotProduct(horizPointing);
		
		realIntegrand[n] *= (beam * cos(2 * pi * freq * 1e6 / c * baselineDotkHat));
		imagIntegrand[n] *= (beam * sin(2 * pi * freq * 1e6 / c * baselineDotkHat));
	}

	//Integrate over the whole sphere and retrun the result;
	xcomplex<double> visibility(0.0,0.0);
	for (int n = 0; n < N; n++){
		visibility.real() += realIntegrand[n];
		visibility.imag() += imagIntegrand[n];
		if (isnan(visibility.real()) != 0 || isnan(visibility.imag()) != 0){
			cout << endl << "NaN Visibility Detected at Pixel " << n << endl << endl;
			break;	
		}
	}

	return visibility * 4 * pi / N;	
}

//Interpolate between calculated visibilities to get visibilities at the desired times
vector<xcomplex<double> > interpolateVisibilities(vector<double>& LSTtimesInterpolated, vector<double>& LSTtimes, vector<xcomplex<double> >& visibilities){
	vector<xcomplex<double> > interpolatedVisibilities;
	int counterOrig = 0;
	for (int n = 0; n < LSTtimesInterpolated.size(); n++){
		while (LSTtimes[counterOrig+1] < LSTtimesInterpolated[n] && counterOrig < (LSTtimes.size()-1)) counterOrig++;
		interpolatedVisibilities.push_back(visibilities[counterOrig] + (visibilities[counterOrig+1] - visibilities[counterOrig]) / (LSTtimes[counterOrig+1] - LSTtimes[counterOrig]) * (LSTtimesInterpolated[n] - LSTtimes[counterOrig]));
	}
	return interpolatedVisibilities;
}

//Saves the results to file with a very detailed filename
void saveResultsToFile(vector<double> LSTtimes, vector<xcomplex<double> > visibilities, int res){
	ofstream outFile;
	outFile.precision(16);
	stringstream ss;
	//ss << resultsFolder << "Revised_Location_Visibilties_for_" << south << "_m_south_" << east << "_m_east_" << up << "_m_up_" << polarization << "_pol_" << freq << "_MHz.dat";
	ss << resultsFolder << "Revised_Location_Visibilties_for_" << round(south) << "_m_south_" << round(east) << "_m_east_" << round(up) << "_m_up_" << polarization << "_pol_" << freq << "_MHz.dat";
	string outFileName = ss.str();
	cout << "Now writing to " << outFileName << endl;
	outFile.open(outFileName.c_str(),ios::trunc);
	for (int n = 0; n < LSTtimes.size(); n++){
		outFile << LSTtimes[n] << " " << visibilities[n].real() << " " << visibilities[n].imag() << endl;
	}
	outFile.close();
}

int main(int argc, char *argv[]){ 
	cout << "NOTE: PRIMARY BEAM TAKEN TO BE A GAUSSIAN WITH FWHM = 10 DEGREES!!!" << endl << endl;

	// Load in all parameters
	loadSpecs();
	if (argc == 4){
		south = atof(argv[1]);
		east = atof(argv[2]);
		up = atof(argv[3]);
		startingLST = 0.0;
		duration = 24;
		timeStep = 10.0/60.0;
	} else {
		cout << endl << "Wrong number of arguments detected.  Please use the python wrapper to run this program.";
		return 1;
	} 
	cout << "BEGIN CALCULATION OF BASELINE (" << south << " m South, " << east << " m East, " << up << " m Up)." << endl << endl << endl;

	string Polarizations[] = {"xx", "yy"};
	//double Frequencies[] = {125.195, 125.586, 125.977, 126.367, 126.758, 127.148, 127.539, 127.93, 128.32, 128.711, 129.102, 129.492, 129.883, 130.273, 130.664, 131.055, 131.445, 131.836, 132.227, 132.617, 133.008, 133.398, 133.789, 134.18, 134.57, 134.961, 135.352, 135.742, 136.133, 136.523, 136.914, 137.305, 137.695, 138.086, 138.477, 138.867, 139.258, 139.648, 140.039, 140.43, 140.82, 141.211, 141.602, 141.992, 142.383, 142.773, 143.164, 143.555, 143.945, 144.336, 144.727, 145.117, 145.508, 145.898, 146.289, 146.68, 147.07, 147.461, 147.852, 148.242, 148.633, 149.023, 149.414, 149.805, 150.195, 150.586, 150.977, 151.367, 151.758, 152.148, 152.539, 152.93, 153.32, 153.711, 154.102, 154.492, 154.883, 155.273, 155.664, 156.055, 156.445, 156.836, 157.227, 157.617, 158.008, 158.398, 158.789, 159.18, 159.57, 159.961, 160.352, 160.742, 161.133, 161.523, 161.914, 162.305, 162.695, 163.086, 163.477, 163.867, 164.258, 164.648, 165.039, 165.43, 165.82, 166.211, 166.602, 166.992, 167.383, 167.773, 168.164, 168.555, 168.945, 169.336, 169.727, 170.117, 170.508, 170.898, 171.289, 171.68, 172.07, 172.461, 172.852, 173.242, 173.633, 174.023, 174.414, 174.805, 175.196, 175.587, 175.978, 176.369, 176.76, 177.151, 177.542, 177.933, 178.324, 178.715, 179.106, 179.497, 179.888, 180.279};
	//double Frequencies[] = {175.1953125, 175.5859375, 175.9765625, 176.3671875, 176.7578125, 177.1484375, 177.5390625, 177.9296875, 178.3203125, 178.7109375, 179.1015625, 179.4921875, 179.8828125, 180.2734375};
	//double Frequencies[] = {125.195, 150.195, 175.196};
	double Frequencies[] = {150.195};

	//Loop over all Polarizations and all Frequencies
	for (int p = 1; p < 2; p++){
		for (int f = 0; f < 1; f++){
			polarization = Polarizations[p];
			freq = Frequencies[f];

			//Figure out which GSM to use and then load its principal components, then properly combine and weight them
			double wavenumber = 2.0 * pi * freq*1e6 / c;
			int res = pickResolution(wavenumber);
			//res = 128;
			vector< Healpix_Map<double> > GSMpcomps = loadGSMpcomps(res);
			Healpix_Map<double> GSM = computeGSM(GSMpcomps, res);

			// Load in all primary beam files for xPol and yPol
			cout << "Now loading primary beam data from file..." << endl;
			char polX[] = "x";
			vector<vector<vector<double> > > beamDataX = makeDataVectorFromFandPol(freq, polX);
			char polY[] = "y";
			vector<vector<vector<double> > > beamDataY = makeDataVectorFromFandPol(freq, polY);

			//Calculate visibility at each time step
			cout << "Now computing visibilities..." << endl;
			double deltaT = timeStep/60;
			if (timeStep/60.0 < 24.0/(res*10+1)) deltaT = 24.0/(res*10+1);	
			vector<xcomplex<double > > visibilities;
			vector<double> LSTtimes;
			for (double t = startingLST; t < (startingLST + duration + deltaT); t += deltaT){
				LST = t;
				LSTtimes.push_back(LST);
				xcomplex<double> visibility= computeVisibility(GSM,beamDataX,beamDataY,res);		
				visibilities.push_back(visibility);
				cout << " " << floor((LST - startingLST)/duration*100) << "% done. \r" << std::flush;
			}
			
			//Interpolate to get the result at the desired timesteps
			if (timeStep/60 < 24.0/(res*10+1)){
				cout << "Now interpolating to get the desired times..." << endl;
				vector<double> LSTtimesInterpolated;
				for (double t = startingLST; t < (startingLST + duration); t += timeStep/60.0){
					LSTtimesInterpolated.push_back(t);
				}
				vector<xcomplex<double> > visibilitiesInterpolated = interpolateVisibilities(LSTtimesInterpolated, LSTtimes, visibilities);
				visibilities = visibilitiesInterpolated;
				LSTtimes = LSTtimesInterpolated;
			}

			//Save the results and termiante
			cout << endl;
			saveResultsToFile(LSTtimes,visibilities,res);
		}
	}

	cout << "Program executed successfully." << endl << endl;
	return 0;
}

