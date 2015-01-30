import ephem
import os
import sys

#Load in specifications of this visibility calculation
#ARGUEMENTS: south, east, up, freq, startTime, startDate, durationInHours
if (len(sys.argv) != 10):
	print '\nToo few or too many arguments detected. The default format is:'
	print 'b_south, b_east, b_up, pol, freq, YYYY/MM/DD, HH:MM:SS, duration, timeStep (in m, m, m, (xx/xy/yy), MHz, date/time UTC, hours, minutes)'
	south = 6
	east = 6
	up = 0
	pol = 'xx'
	freq = 125
	dateAndTime = '2013/07/04 00:00:00'
	UTCtime = '00:00:00'
	duration = 1
	timeStep = 5
	print 'Default set to (6 m, 6 m, 0 m), xx, 125 MHz, 2013/07/04 00:00:00, 1 hour, 5 minutes\n'
else:
	south = float(sys.argv[1])
	east = float(sys.argv[2])
	up = float(sys.argv[3])
	pol = sys.argv[4]
	freq = float(sys.argv[5])
	dateAndTime = sys.argv[6] + ' ' + sys.argv[7]
	UTCtime = sys.argv[7]
	duration = float(sys.argv[8])
	timeStep = float(sys.argv[9])

splitUTCTime = UTCtime.split(':')
UTCTimeDecimal = float(splitUTCTime[0]) + float(splitUTCTime[1])/60 + float(splitUTCTime[2])/60/60

#Load in specifications from the specs file
with open('Specifications.txt') as f:
    specs = f.readlines()

#Calculate local sidereal time
ArraySite = ephem.Observer();
ArraySite.lat = specs[16]
ArraySite.long = specs[19]
ArraySite.date = dateAndTime

splitSideTime = str(ArraySite.sidereal_time()).split(':')
sideTime = float(splitSideTime[0]) + float(splitSideTime[1])/60 + float(splitSideTime[2])/60/60

os.system('make')
dateAndTimeWithPeriods = (dateAndTime.replace("/",".")).replace(":",".")
cmd = './GSMVisibilities' + ' ' + str(south) + ' ' + str(east) + ' ' + str(up) + ' ' + pol + ' ' + str(freq) + ' ' + str(sideTime) + ' ' + str(duration) + ' ' + str(timeStep) + ' ' + str(UTCTimeDecimal) + ' ' + dateAndTimeWithPeriods
print cmd
os.system(cmd)

