"""
File: 'script.py'
The preliminary workings for a script that will take an EPOP RRI HDF5 file and
look at its geographic location in the data file, infer the nearby* SuperDARN
radars for the duration of the dataset, and glean their setting data while 
tagging or making note of the relevant SuperDARN data sets themselves.
Author: David Fairbairn
Date: May 2016
"""

import os
import subprocess

import logging

import davitpy
from davitpy import pydarn
from datetime import datetime
import numpy as np

import timeit
import math

"""
The user can provide an argument specifying the RRI data file, or by default
the script will use the RRI file for April 1st of 2016.

If an argument is provided alongside the script specifying an RRI file, 

"""
import sys
if sys.argv.__len__() == 2 and isinstance(sys.argv[1], str):
    dat_fname = sys.argv[1]
else:
    print "No RRI file specified - going with default..."
    dat_fname = "/home/david/pyth_stuff/script/data/RRI_20160401_072714_073111_lv1_v1.h5" # An RRI data file

# Remotely accessing data on Maxwell.
# The sshfs tool is used to mount Maxwell's data directory locally.
# First, the script unmounts anything currently already mounted in the mounting
# directory (using os.system so that if nothing is there it doesnt crash the 
# script).

"""
import os
os.system("fusermount -uq ./data/remote/")
print "Accessing data files on maxwell, enter your password: "
output = subprocess.check_output(["sshfs", "fairbairn@maxwell.usask.ca:/data/","./data/remote"])
"""

# Some extra code tidbits that can be deleted when I want

# Creating an FOV
site = pydarn.radar.site(code='inv')
myFov = pydarn.radar.radFov.fov(site=site,altitude=300.0,model='IS',coords='geo',ngates=75)

# Taking Long/Lat values from corners of the FOV
rlons,rlats=(np.array(myFov.lonFull)+360.)%360.0,np.array(myFov.latFull)
# np.shape((rlons,rlats)) #(2,17,76) 


"""
---------------------------------- PART 1 --------------------------------------
                        PARSING HDF5 EPHEMERIS DATA
--------------------------------------------------------------------------------
First, the script takes the HDF5 file of interest and grabs its Ephemeris data.


"""

# Running commands to grab Ephemeris Data from RRI HDF5 file.

# Extracting longitude and latitude data from RRI ephemeris.
import h5py
f = h5py.File(dat_fname)
geog_longs = f['CASSIOPE Ephemeris']['Geographic Longitude (deg)'].value
geog_lats  = f['CASSIOPE Ephemeris']['Geographic Latitude (deg)'].value
ephem_times = f['CASSIOPE Ephemeris']['Ephemeris MET (seconds since May 24, 1968)'].value

from script_utils import *
times = ephems_to_datetime(ephem_times)

print "First Geographic Longitude: " + str(geog_longs[0])
print "First Geographic Latitude: " + str(geog_lats[0])


"""
---------------------------------- PART 2 --------------------------------------
                        DETERMINING RELEVANT RADARS 
--------------------------------------------------------------------------------
Next, rough math is performed to check if each latitude point corresponds with
any of the radar stations.

#operations: #radar sites x #lat/lon points in RRI data
"""



# Approach 1: brute force by running Angeline's code on each operational radar 
# for the given latitude-longitude pair
import rgCoords

""" # Using script to go through 11 lat/lon pairs with Angeline Code
# Now, spin through each operational radar and test if it reaches (0,0)
start = timeit.default_timer() # For showing timing
nw = pydarn.radar.network()
results = dict()
lat_subset = geog_lats[0:10]
lon_subset = geog_longs[0:10]
for rad in nw.radars:
    if rad.status == 1:
        bm,gt,view=rgCoords.pos_to_rg(lat_subset,lon_subset,alt=300.0,
coords="geo",dtime=datetime(2016,04,01),frang=180.0,rsep=45.0,rad_id=rad.id) 
        results[rad.name]=(bm,gt,view)

for item in results.items():
    print item

# Takes 49 seconds doing the full brute-force calculations with input of (0,0)
# show results. For (0,0) input, no radar should reach it (maybe azores could).
# Doing 11 samples from Ephemeris data, takes 65.7 seconds.
"""

# Approach 2: brute force by running Ashton's code twice per radar (for both 
# the front and back views)
import range_cells
nw = pydarn.radar.network()
results = dict()
lat_subset = geog_lats#[0:10]
lon_subset = geog_longs#[0:10]
#lat_subset = [0]
#lon_subset = [0]
relevant_radars = dict()

start_t = timeit.default_timer()
for rad in nw.radars:
    # Show progress on screen
    sys.stdout.write(".")
    sys.stdout.flush()

    if rad.status == 1:
        fov_f = davitpy.pydarn.radar.radFov.fov(site=rad.sites[0],altitude=300.,
            model='IS',coords='geo',ngates=75,fov_dir='front')
        bm_f,gt_f=range_cells.findRangeCell(lat_subset,lon_subset,fov_f)        

        fov_b = davitpy.pydarn.radar.radFov.fov(site=rad.sites[0],altitude=300.,
            model='IS',coords='geo',ngates=75,fov_dir='back')
        bm_b,gt_b= range_cells.findRangeCell(lat_subset,lon_subset,fov_b)

        results[rad.name]=(bm_f,gt_f,bm_b,gt_b)

        # non_nan_X contains indices of items that aren't nan (ie. within fov)
        non_nan_f = [n for n in range(np.size(bm_f)) if not math.isnan(bm_f[n])]        
        non_nan_b = [n for n in range(np.size(bm_b)) if not math.isnan(bm_b[n])]
        
        if (non_nan_b.__len__() > 0):
            start = non_nan_b[0] 
            end = non_nan_b[ non_nan_b.__len__() - 1 ] 
            relevant_radars[rad.name] = (rad.code[0],"back", bm_b[start],gt_b[start],bm_b[end],gt_b[end])

        if (non_nan_f.__len__() > 0):
            start = non_nan_f[0]
            end = non_nan_f[ non_nan_f.__len__() - 1 ]
            relevant_radars[rad.name] = (rad.code[0],"front",bm_f[start],gt_f[start],bm_f[end],gt_f[end])

end_t = timeit.default_timer()
print "\nTime required to compute detailed intersections by brute force: " + \
    str(end_t - start_t) + " seconds."

# In general there will be a post-processing step in which only the conjunction
# results will be included in the output

# Output results to show things off:
print "Start time: " + str(ephem_to_datetime(ephem_times[0]))
print "End time: " + str(ephem_to_datetime(ephem_times[-1]))
print "Radar | Fov (f or b), beam_f, gate_f, beam_b, gate_b"
for r in relevant_radars:
    print r + ': ' + str(relevant_radars[r])



"""
---------------------------------- PART 3 --------------------------------------
                      GRABBING RELEVANT SUPERDARN DATA             
--------------------------------------------------------------------------------
Having determined radars, beams, gates, grab the relevant data files (errlog &
perhaps more?) to include here.

Generate or grab a plot showing the satellite's track overlaid vs the superdarn
beams?  

"""


# To convert Ephemeris time data (which is measured since May 24, 1968, we must
# count the time difference between May 24 1968 and 'Epoch' date (Jan 1st 1970)
# which we do by using the "date" command in Bash.
#
# Then, after calculating this offset, we can add the offset to the Ephemeris
# times and convert to UTC using a function from the datetime module.
#
# At first the script accomplished this manually, but now I have put these
# commands into their own function from script_utils, 'ephem_to_datetime'.

from script_utils import *
"""
t2 = subprocess.check_output(["date", "--date=1968-05-24 0:00:00", "+%s"])
t_epoch_eph = int(t2.split("\n",1)[0])

t_sec = t_epoch_eph + int(ephem_times[0]) # first ephemeris timing element
start = datetime.utcfromtimestamp(t_sec)

t_sec = t_epoch_eph + int(ephem_times[-1]) # last ephemeris timing element
end = datetime.utcfromtimestamp(t_sec)
"""
start = ephem_to_datetime(ephem_times[0])
end = ephem_to_datetime(ephem_times[-1])


stmonth = "0" + str(start.month) if str(start.month).__len__() == 1 else str(start.month)
stday = "0" + str(start.day) if str(start.day).__len__() == 1 else str(start.day)

# TODO: automated selection of radar codes to replace 'rkn' below
fname = str(start.year) + str(stmonth) + str(stday) + "." + "rkn.errlog.bz2"
tstring = str(start.hour) + ":" +  str(start.minute) + ":00"


# The SuperDARN errlog files are compressed in .bz2 file formats. However, the
# bz2 library provides functions for reading .bz2 files like normal text files.
"""
# TODO:automated selection of radar codes to replace 'rkn' below
import  bz2
fdir = "rkn" + "_errlog"
f = bz2.BZ2File("./data/remote/rkn_errlog/" + fname)

# With the file open, search for the desired time interval's beginning.
# #TODO: (and possibly its ending)

found = False
while not found:
    ln = f.readline()
    if ln.find(tstring) != -1:
        found = True
        num = f.tell()
        print str(num) + ": " + ln

# Having determined the relevant line of text in the errlog file, grab all the
# relevant errlog data spanning the ephemeris file
f.seek(num) 
rel_lines = f.readline()
# TODO: explicitly determine the required duration. (search for end time)

for i in range(500): # 500 should cover the 237 data samples and in general, a 4 minute pass
    rel_lines = rel_lines + f.readline()
"""
uofs_rads = []
for r in relevant_radars:
    if r in ['Saskatoon','Prince George','Clyde River','Inuvik','Rankin Inlet']:
        uofs_rads.append(r) 
        plot_fov_sat(r,start,geog_longs,geog_lats)
