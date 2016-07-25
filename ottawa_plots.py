
import matplotlib.pyplot as plt
from script_utils import *

import math
import numpy as np
import sys



# **** CHOOSE ONE OF THESE RRI FILES THEN RUN THE SCRIPT ****
geog_longs,geog_lats,ephemtimes = get_rri_ephemeris("./data/RRI_20160418_222759_223156_lv1_v2.h5") #18th
index_inversion = 167 #for 18th

#index_inversion = 178 #for 19th
#geog_longs,geog_lats,ephemtimes = get_rri_ephemeris("./data/RRI_20160419_220939_221336_lv1_v2.h5") #19th

#index_inversion = 213 #for 20th
#geog_longs,geog_lats,ephemtimes = get_rri_ephemeris("./data/RRI_20160420_215117_215514_lv1_v2.h5") #20th


ottawa_long = (-75.6972+360.)%360.0 # Make it positive valued
ottawa_lat = 45.4215
times = ephems_to_datetime(ephemtimes)

geog_longs = (geog_longs+360.)%360.0

dists = []
longdists = []
latdists = []
shortest_dist = sys.maxint
for i in range(np.size(geog_longs)):
    dist = haversine(geog_longs[i], geog_lats[i], ottawa_long, ottawa_lat)
    longdist = abs(geog_longs[i] - ottawa_long)
    latdist = abs(geog_lats[i] - ottawa_lat)
    #dist = np.sqrt(longdist*longdist + latdist*latdist)  
    if dist < shortest_dist:
        shortest_dist = dist
        shortest = i
    dists.append(dist)
    longdists.append(longdist)
    latdists.append(latdist)

appr_time = times[shortest]

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
