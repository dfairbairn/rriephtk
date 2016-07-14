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
import sys

import logging #TODO: use this.

import davitpy
from davitpy import pydarn
from datetime import datetime
import numpy as np

import timeit
import math


# -------------------------------- PART 1 --------------------------------------
#                              INITIALIZATION 
# ------------------------------------------------------------------------------
# First, the script determines which RRI data it's looking at.

# The user can provide an argument specifying the RRI data file, or by default
# the script will use the RRI file for April 2nd of 2015.
#

# Test if directory structure with ./data, ./data/output, ./data/remote exist
if not os.path.isdir('./data'):
    os.system('mkdir ./data')
    print "Creating script directory 'data/'."
if not os.path.isdir('./data/output'):
    os.system('mkdir ./data/output')
    print "Creating script directory 'data/output/'."
if not os.path.isdir('./data/remote'):
    os.system('mkdir ./data/remote')
    print "Creating script directory 'data/remote/'."

import sys
if sys.argv.__len__() == 2 and isinstance(sys.argv[1], str):
    dat_fname = sys.argv[1]
else:
    print "No RRI file specified - going with default..."
    dat_fname = "./data/RRI_20150402_032244_033241_lv1_v2.h5" # An RRI data file
    
if not os.path.exists(dat_fname):
    print "No RRI file by that name. Exitting."
    exit()

# The sshfs tool is used to mount Maxwell's data directory locally.
# 
# First, the script unmounts anything currently already mounted in the mounting
# directory (using os.system so that if nothing is there it doesnt crash the 
# script).
os.system("fusermount -uq ./data/remote/")
print "Accessing data files on maxwell, enter your username: "
username = raw_input("Username: ")
print "Accessing data files on maxwell, enter your password: "
output = subprocess.check_output(["sshfs", username + "@maxwell.usask.ca:/data/","./data/remote"])


# Some extra code tidbits that can be deleted when I want

# Creating an FOV
site = pydarn.radar.site(code='inv')
myFov = pydarn.radar.radFov.fov(site=site,altitude=300.0,model='IS',coords='geo',ngates=75)

# Taking Long/Lat values from corners of the FOV
rlons,rlats=(np.array(myFov.lonFull)+360.)%360.0,np.array(myFov.latFull)
# np.shape((rlons,rlats)) #(2,17,76) 



# -------------------------------- PART 2 --------------------------------------
#                       PARSING HDF5 EPHEMERIS DATA
# ------------------------------------------------------------------------------
# First, the script takes the HDF5 file of interest and grabs its Ephemeris data.

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


# -------------------------------- PART 3 --------------------------------------
#                        DETERMINING RELEVANT RADARS 
# ------------------------------------------------------------------------------
# Next, rough math is performed to check if each latitude point corresponds with
# any of the radar stations.

# Current approach: brute force by running Ashton's code twice per radar (for 
# both the front and back views). Used to use Angeline Burrell's code.
import range_cells
nw = pydarn.radar.network()
results = dict()
lat_subset = geog_lats#[0:10]
lon_subset = geog_longs#[0:10]
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



# -------------------------------- PART 4 --------------------------------------
#                 GRAB SUPERDARN DATA, CREATE OUTPUT FILES             
# ------------------------------------------------------------------------------
# Having determined radars, beams, gates, grab the relevant data files (errlog &
# perhaps more?) to include here.



# To convert Ephemeris time data (which is measured since May 24, 1968, we must
# count the time difference between May 24 1968 and 'Epoch' date, Jan 1st 1970)
# which we do by using the "date" command in Bash.
#
# Then, after calculating this offset, we can add the offset to the Ephemeris
# times and convert to UTC using a function from the datetime module.
#
# At first the script accomplished this manually, but now I have put these
# commands into their own function from script_utils, 'ephem_to_datetime'.

from script_utils import *

start = ephem_to_datetime(ephem_times[0])
end = ephem_to_datetime(ephem_times[-1])

# Check if they're radars for which we actually have SuperDARN data currently 
uofs_rads = []
for r in relevant_radars:
    if r in ['Saskatoon','Prince George','Clyde River','Inuvik','Rankin Inlet']:
        uofs_rads.append(r) 

# ERRLOG files contain entries describing the beam number and frequency of 
# the transmission, occurring roughly every three seconds, starting on each 
# minute (at the ":00" mark).
#
# FURTHERMORE: the fact that the timestamps on SuperDARN data are untrustworthy
# means we need to look at SuperDARN information up to a few minutes before and
# after the RRI times of interest.
#
# Thus, this script will return all ERRLOG data starting *3* minutes before
# the RRI's first timestamp, and ending *3* minutes after the RRI's last
# timestamp, so as to give us some margin of error. The frequency that were
# transmitted upon will hopefully inform the user as to where the two time
# logs actually sync up.


# Also, need to pad the datetime objects so that their format matches the 
# filename and timestamp formats for the ERRLOG files.
st_month = "0" + str(start.month) if str(start.month).__len__() == 1 else str(start.month)
st_day = "0" + str(start.day) if str(start.day).__len__() == 1 else str(start.day)
st_hour = "0" + str(start.hour) if str(start.hour).__len__() == 1 else str(start.hour)
end_hour = "0" + str(end.hour) if str(end.hour).__len__() == 1 else str(end.hour)
# Also, due to untrustworthy timestamps, we look at times in the ERRLOG file
# up to 3 minutes before and after the RRI times.
st_min = start.minute-3
end_min = end.minute+3
st_min = "0" + str(st_min) if str(st_min).__len__() == 1 else str(st_min)
end_min = "0" + str(end_min) if str(end_min).__len__() == 1 else str(end_min)


start_string = str(st_hour) + ":" +  str(st_min) + ":" # "00" # Omit seconds for now in case
end_string = str(end_hour) + ":" + str(end_min) + ":" #"00" # they don't follow expected pattern
# NOTE: THIS DOESNT WORK FOR RKN_ERRLOG FOR 20150402 BECAUSE END_STRING DOESNT APPEAR
# TODO: MAKE THE ENDPOINT FOR THE ERRLOG READ JUST BE WHEN THE EPOPSOUND PROGRAM CEASES??
#       OR PERHAPS JUST LOOK FOR A TIME 'GREATER' THAN END_TSTRING?

# The SuperDARN errlog files are compressed in .bz2 file formats. However, the
# bz2 library provides functions for reading .bz2 files like normal text files.
import  bz2

for u in uofs_rads:
    rcode = (nw.getRadarByName(u)).code[0]
    plot_fov_sat(u,start,geog_longs,geog_lats,suppress_show=True) 

    fname_bz2 = str(start.year) + str(st_month) + str(st_day) + "." + rcode + ".errlog.bz2"
    fname_reg = str(start.year) + str(st_month) + str(st_day) + "." + rcode + ".errlog"
    fpath = "./data/remote/" + rcode + "_errlog/"

    if os.path.exists(fpath + fname_bz2):
        f = bz2.BZ2File(fpath + fname_bz2)
        logging.debug("Uses bz2 errlog file.")
    elif os.path.exists(fpath + fname_reg):
        logging.debug("Uses non-bz2 errlog file.")
        f = open(fpath + fname_reg)
    else:
        logging.error("No ERRLOG file found!")
        exit()

    # With the file open, search for the desired time interval's beginning.
    # Search for the position of the first line of interest
    found = False
    while not found:
        ln = f.readline()
        if ln.find(start_string) != -1:
            found = True
            start_line = f.tell()
            print str(start_line) + ": " + ln
    # Now search for the position of the final line of interest
    found = False
    while not found:
        ln = f.readline()
        if ln.find(end_string) != -1:
            found = True
            end_line = f.tell()
            print str(end_line) + ": " + ln
        elif f.tell > (start_line + 1000000): #1000000 bytes
            end_line = start_line + 1000000
            found = True
            print str(end_line) + ": " + ln

    # Having determined the relevant line of text in the errlog file, grab all the
    # relevant errlog data spanning the ephemeris file
    f.seek(start_line) 
    rel_lines = f.readline()
    while f.tell() <= end_line:
        rel_lines = rel_lines + f.readline()

    outp = open("./data/output/"+(str(start.year)+st_month+st_day+"_"+st_hour+st_min+"_"+rcode+".dat"),"w+")
    outp.write("OUTPUT FILE FOR RRI CONJUNCTION SCRIPT\n======================================\n")
    outp.write("Start Time: "+start.__str__()+"\nEnd Time: "+end.__str__())
    outp.write("\nRRI File: "+dat_fname+"\n")
    
    for r in relevant_radars:
        outp.write("\n"+str(r) + ": " + str(relevant_radars[r]))
    outp.write("\nFront FOV and Back FOV (beam,gate) pairs for ephemeris position:\n")
    outp.write((results[u]).__str__())
    outp.write("\n\nRelevant SuperDARN data from +/- 3 minutes of the RRI data:\n")
    outp.write(rel_lines)

os.system("fusermount -uq ./data/remote/")

