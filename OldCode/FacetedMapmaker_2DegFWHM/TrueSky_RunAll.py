import os
import sys
import math
import fnmatch

allTests = ("Snapshot_Time/", "PSF_Size/", "Spatial_Resolution/", "Facet_Size/", "Very_Small_Facet/", "Only_Zenith/");

testDirectoryPrefix = "./Test_"
for test in allTests:
	testDirectory = testDirectoryPrefix + test
	for file in os.listdir(testDirectory):
	    if fnmatch.fnmatch(file, 'Spec*'):
	        #cmd = "qsub -cwd -v SPECSFILE='" + testDirectory + str(file) + "',LOGFILE='" + testDirectory + logname + "' -l h_vmem=4G ./singleGridEngineCore.sh"
	       	cmd = "./True_Faceted_Sky/True_Faceted_Sky " + testDirectory + str(file) 
	       	print cmd
	       	os.system(cmd)


