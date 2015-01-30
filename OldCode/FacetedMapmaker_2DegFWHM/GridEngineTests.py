import os
import sys
import math
import fnmatch

allTests = ("Snapshot_Time/", "PSF_Size/", "Spatial_Resolution/", "Facet_Size/", "Very_Small_Facet/", "Only_Zenith/");

resultsDirectoryPrefix = "/nfs/pritchard/r1/jsdillon/Mapmaking_Results/2Deg_FWHM/Results_"
testDirectoryPrefix = "./Test_"
for test in allTests:
	resultsDirectory = resultsDirectoryPrefix + test
	testDirectory = testDirectoryPrefix + test
	for file in os.listdir(testDirectory):
	    if fnmatch.fnmatch(file, 'Spec*'):
	        logname = "log_" + str(file)
	        cmd = "qsub -cwd -v SPECSFILE='" + testDirectory + str(file) + "',LOGFILE='" + resultsDirectory + logname + "' -l h_vmem=4G ./singleGridEngineCore.sh"
	       	print cmd
	       	os.system(cmd)
	       	cpcmd = "cp " + testDirectory + str(file) + " " + resultsDirectory
	       	os.system(cpcmd)