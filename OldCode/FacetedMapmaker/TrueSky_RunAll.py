import os
import sys
import math
import fnmatch

#test = "Snapshot_Time/"
#test = "PSF_Size/"
#test = "Spatial_Resolution/"
#test = "Facet_Size/"
#test = "Very_Small_Facet/"
test = "Only_Zenith/"

testDirectoryPrefix = "./Test_"
testDirectory = testDirectoryPrefix + test

for file in os.listdir(testDirectory):
    if fnmatch.fnmatch(file, 'Spec*'):
        #cmd = "qsub -cwd -v SPECSFILE='" + testDirectory + str(file) + "',LOGFILE='" + testDirectory + logname + "' -l h_vmem=4G ./singleGridEngineCore.sh"
       	cmd = "./True_Faceted_Sky/True_Faceted_Sky " + testDirectory + str(file) 
       	print cmd
       	os.system(cmd)


