import os
import sys
import time

#To call: nohup python -u OmniscopeSimulation.py > & ./Logs/Python_Log.txt &

#Load in antenna positions from the array file
uniqueBaselines = [i.strip().split() for i in open("IntendedUniqueBaselines.dat").readlines()]

for b in range(0, len(uniqueBaselines)):
	#cmd = 'nohup nice -n 15 ./GSMVisibilitiesAllFreqAllPol ' + str(uniqueBaselines[b][0]) + ' ' + str(uniqueBaselines[b][1]) + ' 0.0 >& ./Logs/log_' + str(b) + '.txt &'
	cmd = 'nohup nice -n 15 ./GSMVisibilitiesAllFreqAllPol ' + str(uniqueBaselines[b][0]) + ' ' + str(uniqueBaselines[b][1]) + ' ' + str(uniqueBaselines[b][2]) + ' >& ./Logs/log_' + str(b) + '.txt &'
	print cmd
	os.system(cmd)
	time.sleep(5)
	process = os.popen('ps -C GSMVisibilitiesAllFreqAllPol | grep -c GSMVis')
	nProcRunning = int(process.read())
	while (nProcRunning >= 8):		
		process = os.popen('ps -C GSMVisibilitiesAllFreqAllPol | grep -c GSMVis')
		nProcRunning = int(process.read())
		print 'Waiting to start baseline '+str(b)+'...'
		time.sleep(60)
