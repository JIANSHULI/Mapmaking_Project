targetDir = './Test_Facet_Size/';

facetSize = 12;
angularResolution = .5;
centerPixRA = 266.2417;
centerPixDec = -30.706267;
cornerPixRA = centerPixRA - facetSize/2;
cornerPixDec = centerPixDec - facetSize/2;

raList = {};
decList = {};

for nSide = 1:4
    raList{nSide} = [];
    decList{nSide} = [];
    for a = 1:nSide
        for d = 1:nSide
            filename = [targetDir 'Specs_' num2str(nSide) '_' num2str(a) '_' num2str(d) '.txt'];
            fileID = fopen(filename,'w');
            fprintf(fileID,'frequency = 150.195\n');
            fprintf(fileID,'polarization = yy\n');
            fprintf(fileID,'baselinesFile = ../GSM_Visibilities/IntendedUniqueBaselines.dat\n');
            fprintf(fileID,'visibilitiesFolderAndPrefix = ../GSM_Visibilities/Revised_Location_Visibilties_for_\n');
            fprintf(fileID,'unnormalizedFullMapFilename = unnormalized_full_map.dat\n');
            fprintf(fileID,'finalMapFilename = map.dat\n');
            fprintf(fileID,'PSFfilename = PSF.dat\n');
            fprintf(fileID,'noiseCovFilename = noiseCov.dat\n');
            fprintf(fileID,'DmatrixFilename = Dmatrix.dat\n');
            fprintf(fileID,'dataProductFilePrefix = /nfs/pritchard/r1/jsdillon/Mapmaking_Results/Results_Facet_Size/%d_%d_%d_\n',nSide,a,d);
            fprintf(fileID,'arrayLatitude = -30.706267\n');
            fprintf(fileID,'arrayLongitude = 21.321317\n');
            ra = cornerPixRA + facetSize/nSide/2 + (a-1)*facetSize/nSide;
            dec = cornerPixDec + facetSize/nSide/2 + (d-1)*facetSize/nSide;
            raList{nSide}(a,d) = ra;
            decList{nSide}(a,d) = dec;
            fprintf(fileID,'facetRA = %f\n',ra);
            fprintf(fileID,'facetDec = %f\n',dec);
            fprintf(fileID,'facetSize = %f\n',facetSize/nSide);
            fprintf(fileID,'angularResolution = %f\n',angularResolution);
            psfFactor = ceil(3*nSide/2)*2+1; 
            fprintf(fileID,'PSFextensionBeyondFacetFactor = %d\n',psfFactor);
            fprintf(fileID,'snapshotTime = 120\n');
            fprintf(fileID,'integrationTime = 10.0\n');
            fprintf(fileID,'channelBandwidth = 0.048828125\n');
            fprintf(fileID,'gaussianPrimaryBeam = 1\n');
            fprintf(fileID,'primaryBeamFWHM = 10.0\n');
            fprintf(fileID,'xpolOrientationDegreesEastofNorth = 90.0\n');
            fprintf(fileID,'maximumAllowedAngleFromPBCenterToFacetCenter = 10\n');
            fprintf(fileID,'noiseStd = 100\n');
            fclose(fileID);
        end
    end
end
            
