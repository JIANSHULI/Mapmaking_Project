import os
import sys
import math
import fnmatch

test = "MultiFrequency/"


nGigsRAM = 1

resultsDirectoryPrefix = "/nfs/pritchard/r1/jsdillon/Mapmaking_Results/Results_"
testDirectoryPrefix = "./Test_"
resultsDirectory = resultsDirectoryPrefix + test
testDirectory = testDirectoryPrefix + test
for file in os.listdir(testDirectory):
    if fnmatch.fnmatch(file, 'Spec*'):
        logname = "log_" + str(file)
        cmd = "qsub -cwd -v SPECSFILE='" + testDirectory + str(file) + "',LOGFILE='" + resultsDirectory + logname + "' -l h_rt=12:00:00,h_vmem=" + str(nGigsRAM) + "G -pe chost 1 ./singleGridEngineCore.sh"
       	
       	print cmd
       	os.system(cmd)
       	cpcmd = "cp " + testDirectory + str(file) + " " + resultsDirectory
       	os.system(cpcmd)