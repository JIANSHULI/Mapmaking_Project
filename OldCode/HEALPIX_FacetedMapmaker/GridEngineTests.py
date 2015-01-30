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

resultsDirectoryPrefix = "/nfs/pritchard/r1/jsdillon/Mapmaking_Results/Results_"
testDirectoryPrefix = "./Test_"
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