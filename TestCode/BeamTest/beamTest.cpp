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

const double pi = 3.1415926535897932384626433832795;
int nAltBeamPoints = 45; //running from 0 to 90 in increments of 90/44
int nAzBeamPoints = 180; //running for 0 to 360 in increments of 360/179
double firstBeamFreq = 110; //MHz
double beamDeltaFreq = 10; //MHz
double lastBeamFreq = 190; //MHz
string beamDirectory = "./mwa_beam_";


int nFreqFiles = (lastBeamFreq - firstBeamFreq)/10 + 1;
double deltaBeamAlt = 90.0 / (nAltBeamPoints-1);
double deltaBeamAz = 360.0 / (nAzBeamPoints-1);
string polarization;
double freq = 150.195; //MHz
double xpolOrientationDegreesEastofNorth = 90.0;

struct pointing {
	double theta, phi;
	pointing(double thetaIn, double phiIn) : theta(thetaIn), phi(phiIn) {}
};

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
	  cout << numFreqFiles << " freqs.line# " << lineNum<<endl; 
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
	    cout << "Offset is " << offsetFreq<<endl;
	  }
	  else{ // if interpolating gain between frequencies:
	    dataOutN=2;
	    // set offset to be lower frequency block
	    // and copy the next frequency block as well
	    for(n=11;n<20;n++){
	      if((f-(n*10))<10){
	    	offsetFreq=(n-11);//*dataBlockSize;
	    	cout << "Offset for doubleBlock is " << offsetFreq<<endl;
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

//Returns the value of the primary beam for a given pointing in horizontal coordiantes
double primaryBeamShana(pointing& horizPointing, vector<vector<vector<double> > >& beamDataX, vector<vector<vector<double> > >& beamDataY){
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


//This funciton loads the primary beam into a array of discretized values of alt and az.
//The data file, somewhat inconsistently, has azimuth 0 in the direction pointed by the XX polarization and increases CCW
vector< vector<double> > loadDiscretizedPrimaryBeam(){
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

//This function retruns the value of the gain of the PB as a funciton of altitude and azimuth
//Unfortunately, the primary beam azimuth is stored with the XX polarization as azimuth zero and continues CCW. 
double primaryBeam(horizPoint& pointing, vector< vector<double> >& PB){
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

//	double firstBeamFreq = 110; //MHz
//double beamDeltaFreq = 10; //MHz
//double lastBeamFreq = 190; //MHz

int main() {
	cout << endl << "Running Beam Test..." << endl << endl;
	char polX[] = "x";
	vector<vector<vector<double> > > beamDataX = makeDataVectorFromFandPol(freq, polX);
	char polY[] = "y";
	vector<vector<vector<double> > > beamDataY = makeDataVectorFromFandPol(freq, polY);
	
	polarization = "xx";
	
	pointing horizPointinShana(pi/2, pi/8);
	cout << primaryBeamShana(horizPointinShana, beamDataX, beamDataY) << endl;
	//cout << "Done with Shana's version." << endl << endl;
	

	vector< vector<double> > primaryBeamDiscretized = loadDiscretizedPrimaryBeam();
	horizPoint myPointing(pi/2, pi/8);
	cout << primaryBeam(myPointing, primaryBeamDiscretized) << endl;

	cout << "Comparison: " << endl << endl;

	for (double th = 0; th< 2*pi; th+=.1*pi){
		pointing horizPointinShana2(pi/3, th);
		horizPoint myPointing2(pi/3, th);

		cout << primaryBeamShana(horizPointinShana2, beamDataX, beamDataY) << ", " << primaryBeam(myPointing2,primaryBeamDiscretized) << endl;
	}

	cout << "Done." << endl << endl;
	return 0;
}
