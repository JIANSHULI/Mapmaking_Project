import os
import sys
import math
import fnmatch

specs = ["frequency = 150.195",
"polarization = yy",
"baselinesFile = ../GSM_Visibilities/IntendedUniqueBaselines.dat",
"visibilitiesFolderAndPrefix = ../GSM_Visibilities/Revised_Location_Visibilties_for_",
"trueSkyFilename = trueSky.dat",
"finalMapFilename = map.dat",
"PSFfilename = PSF.dat",
"healpixPixelFilename = healpixPixels.dat",
"extendedHealpixPixelFilename = extendedHealpixPixels.dat",
"pixelCoordinatesFilename = pixelCoordinates.dat",
"extendedPixelCoordiantesFilename = extendedPixelCoordinates.dat",
"noiseCovFilename  = noiseCov.dat",
"DmatrixFilename = Dmatrix.dat",
"dataProductFilePrefix = /nfs/pritchard/r1/jsdillon/Mapmaking_Results/Results_MultiFrequency/",
"arrayLatitude = -30.706267",
"arrayLongitude = 21.321317",
"facetRA = 265.4296875",
"facetDec = -30.706267",
"facetSize = 10",
"NSIDE = 1024",
"PSFextensionBeyondFacetFactor = 3",
"snapshotTime = 10",
"integrationTime = 10.0",
"channelBandwidth = 0.048828125",
"gaussianPrimaryBeam = 1",
"primaryBeamFWHM = 10.0",
"xpolOrientationDegreesEastofNorth = 90.0",
"maximumAllowedAngleFromPBCenterToFacetCenter = 5",
"noiseStd = 100",
"overwriteVisibilitiesWithGSM = 1",
"GSMres = 128",
"overwriteVisibilitiesWithPointSources = 0",
"alsoComputePointSourcePSF = 0",
"pointSourceFile = mwacs_all_b3_140206.dat",
"pointSourceReferenceFreq = 180",
"pointSourceFluxUpperLimit = 1",
"pointSourceFluxLowerLimit = .001",
"pointSourceOutputCatalogFilename = pointSources.dat",
"pointSourcePSFFilename = pointSourcePSF.dat",
"throwOutFluxOutsideExtendedFacet = 1"]



f = 80
while f <= 200:
	filename = "Specifications_" + str(f) + ".txt"
	thisFile = open(filename,'w')
	for n in range(len(specs)):
		if "frequency = " in specs[n]:
			thisFile.write("frequency = " + str(f) + "\n")
		elif "dataProductFilePrefix = " in specs[n]:
			thisFile.write(specs[13] + str(f) + "_" + "\n")
		else:
			thisFile.write(specs[n] + "\n")		
	thisFile.close() 
	f += .5
