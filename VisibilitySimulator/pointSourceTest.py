import ephem
import os
import sys

dateAndTime = '2012/05/27 04:07:49'

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
print "LST: ", sideTime

# Latitude and longitude for galactic coordinates;
#ga = ephem.Galactic('202.5','36.4236')
ga = ephem.Galactic('225.703','10.8069')
eq = ephem.Equatorial(ga)

print "Long: ", ephem.degrees(225.703/360*2*3.1415926)
print "Lat: ", ephem.degrees(10.8069/360*2*3.1415926)

bodyString = "Test,f|M|F7," + str(eq.ra) + "," + str(eq.dec) + ",0,2012"
#print bodyString
testPS = ephem.readdb(bodyString)
testPS.compute(ArraySite)
#print testPS.dec
print "RA: ", ephem.degrees(eq.ra)
print "Dec: ", ephem.degrees(eq.dec)
print "Azimuth: ", testPS.az
print "Altitude: ", testPS.alt



#print eq.ra
#print eq.dec

#print ga.long, ga.lat