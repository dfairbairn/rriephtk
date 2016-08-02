"""
file: 'ottawa_plots.py'
author: David Fairbairn
date: July 2016
description: This file contains code used to produce fly-by plots of CASSIOPE
    vs the Ottawa-based radar that was running an RRI experiment on the dates
    April 18-22nd.

    On each day, faraday rotation behaviour was seen to be occurring, reaching
    a point where it would stop rotating in one direction and start rotating
    the opposite way. This turned out to _NOT_ be a simple case of the distance
    from CASSIOPE to Ottawa reaching a minimum at this same point, so these 
    plots serve to contrast the location of the rotation inversion, and the 
    point of closest approach.

    In addition, magnetic coordinates and IGRF can be included.

"""
from script_utils import *
from data_utils import *

from davitpy.utils import *

import math
#import matplotlib.pyplot as plt
#import numpy as np
import sys

# **** CHOOSE ONE OF THESE RRI FILES THEN RUN THE SCRIPT ****
geog_longs,geog_lats,ephemtimes = get_rri_ephemeris("./data/RRI_20160418_222759_223156_lv1_v2.h5") #18th
index_inversion = 167 #for 18th

#index_inversion = 178 #for 19th
#geog_longs,geog_lats,ephemtimes = get_rri_ephemeris("./data/RRI_20160419_220939_221336_lv1_v2.h5") #19th

#index_inversion = 213 #for 20th
#geog_longs,geog_lats,ephemtimes = get_rri_ephemeris("./data/RRI_20160420_215117_215514_lv1_v2.h5") #20th

#index_inversion = 205 #?? for 21st?
#geog_longs,geog_lats,ephemtimes = get_rri_ephemeris("./data/RRI_20160421_213255_213652_lv1_v2.h5") #21st

#index_inversion = 222 #for 22nd
#geog_longs,geog_lats,ephemtimes = get_rri_ephemeris("./data/RRI_20160422_211435_211832_lv1_v2.h5") #22nd

# Location of Ottawa: I looked it up and am hard-coding it here.
ottawa_long = -75.6972
ottawa_lat = 45.4215
times = ephems_to_datetime(ephemtimes)

# Make all longitudes positive?
#ottawa_long = (ottawa_long+360.)%360.0
#geog_longs = (geog_longs+360.)%360.0

# *** FINDING THE CLOSEST APPROACH ***
# Using the Haversine formula (in a function in script_utils.py), the closest
# approach is determined by brute force.
dists = []
longdists = []
latdists = []
shortest_dist = sys.maxint #Initially set this very high
for i in range(np.size(geog_longs)):
    dist = haversine(geog_longs[i], geog_lats[i], ottawa_long, ottawa_lat)
    # Initially I took a quick and dirty approach to find the point of smallest
    # Euclidean distance in terms of latitudes and longitudes.
    #longdist = abs(geog_longs[i] - ottawa_long) # difference of longitudes
    #latdist = abs(geog_lats[i] - ottawa_lat) # difference of latitudes
    #dist = np.sqrt(longdist*longdist + latdist*latdist)  
    if dist < shortest_dist:
        shortest_dist = dist
        shortest = i

appr_time = times[shortest]

# The numeric data type that I was retrieving from geog_longs, when _NOT_ stored
# in an array, was being rejected by the mapObj() function below. So I convert 
# these numbers to floats explicitly here.
shortest_ephem_long = float(geog_longs[shortest])
shortest_ephem_lat = float(geog_lats[shortest])
inversion_ephem_long = float(geog_longs[index_inversion])
inversion_ephem_lat = float(geog_lats[index_inversion])

"""
# From when we weren't bothering to try to print on a background map of the world.
plt.plot(ottawa_long,ottawa_lat,'ro',label="Ottawa")
plt.plot(geog_longs,geog_lats,label="EPOP ground track")
plt.plot([geog_longs[shortest], ottawa_long],[geog_lats[shortest], ottawa_lat],'g')
plt.plot(geog_longs[shortest],geog_lats[shortest],'bo',label=("Shortest Approach at " + str(appr_time)))

plt.plot(geog_longs[index_inversion],geog_lats[index_inversion],'yo',label=("Inversion of Faraday Rotation"))

plt.xlabel('Geographic Longitude (degrees)')
plt.ylabel('Geographic Latitude (degrees)')
plt.title("EPOP Closest Approach vs. Ottawa radar for " + "2016-04-" + str(times[0].day))
plt.legend()
plt.show()
"""
# A different font for the legend etc. might be nice
font = {'fontname':'Computer Modern'}
m = plotUtils.mapObj(lat_0=38.0, lon_0=-76.0, width=111e3*20, height=111e3*50, coords='geo',datetime=times[0])

# FIRST: Plot the location of Ottawa
x,y = m(ottawa_long,ottawa_lat,coords='geo')
m.plot(x,y,'ro',label="Ottawa")

# SECOND: Plot the satellite ground-track.
x,y = m(geog_longs, geog_lats, coords='geo')
m.plot(x,y,'b',label="EPOP ground track")

# THIRD: Plot a circle emphasizing the point of closest approach
x,y = m(shortest_ephem_long,shortest_ephem_lat, coords='geo')
m.plot(x,y,'bo',label=("Shortest Approach at " + str(appr_time)))

# FOURTH: Plot the line from Ottawa to the nearest approach of the satellite.
x,y = m([shortest_ephem_long, ottawa_long], [shortest_ephem_lat, ottawa_lat], coords='geo')
m.plot(x,y,'g')

# FIFTH: Plot the piont I've determined is the point of the Faraday Rotation inversion.
x,y = m(inversion_ephem_long,inversion_ephem_lat,coords='geo')
m.plot(x,y,'yo',label=("Inversion of Faraday Rotation"))

plt.xlabel('Geographic Longitude (degrees)',**font)
plt.ylabel('Geographic Latitude (degrees)',**font)
plt.title("EPOP Closest Approach vs. Ottawa radar for " + "2016-04-" + str(times[0].day),**font)
plt.legend()
plt.show()
