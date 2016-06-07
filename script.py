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

import davitpy
from davitpy import pydarn
from datetime import datetime
import numpy as np

import timeit

#radars = davitpy.pydarn.network()

# Creating an FOV
site = pydarn.radar.site(code='inv')
myFov = pydarn.radar.radFov.fov(site=site,altitude=300.0,model='IS',coords='geo',ngates=75)

# Taking Long/Lat values from corners of the FOV
rlons,rlats=(np.array(myFov.lonFull)+360.)%360.0,np.array(myFov.latFull)
# np.shape((rlons,rlats)) #(2,17,76) 
centerpt = (rlons[0][0],rlats[0][0])
left = (rlons[0][75],rlats[0][75])
right = (rlons[16][75],rlats[16][75])

# Ashton code for extracting beam and gate from lat/lon

# Fixed point using default radar (sas)
#beams,gates = range_cells.findRangeCell([60.0],[-90])
#beams,gates = range_cells.findRangeCell([60.0,64.0],[-90,-100])
#print beams,gates

# Using an array of lats and lons, and for a specific radar ('rkn')
# Note: the lats and lons aren't realistic: they're actually from PyDarn
#import range_cells
#beams,gates = range_cells.findRangeCell(lats,lons, myFov=myFov)


# ANGELINE BURRELL COORDS->RG CODE
import rgCoords

dtime=datetime(2011,06,01)
# rg_gate: range gates, bmnum: beam numbers, fovflg: flag showing forward FOV or back lobe
# NOTE: rlats and rlons currently are actually lat/lon coords meant for 
#       checking FOV containment, not for being contained themselves. JUST A TEST
#rg_gate,bmnum,fovflg = rgCoords.pos_to_rg(rlats[7],rlons[7],coords='geo',alt=300.0,dtime=dtime,frang=180.0,rsep=45.0,rad_code='sas')


"""
---------------------------------- PART 1 --------------------------------------
                        PARSING HDF5 EPHEMERIS DATA
--------------------------------------------------------------------------------
First, the script takes the HDF5 file of interest and grabs its Ephemeris data.


"""

# Running commands to grab Ephemeris Data from RRI HDF5 file.
# 
#

# A now-unnecessary bit of code to generate new data files (and not overwrite)
newfile = "data/output/tmp_ephem.dat" 
tmp_fname = newfile
i = 0
while (os.path.isfile(tmp_fname)):
    i+=1
    tmp_fname = newfile + str(i)

# For now, data file assumed to be in a fixed location. #TODO: custom paths
dat_fname = "/home/david/pyth_stuff/script/data/RRI_20160401_072714_073111_lv1_v1.h5" # An RRI data file

# This bash command is now superfluous, because we use h5py to get data 
bashCmd = "h5dump -d CASSIOPE\ Ephemeris/Geographic\ Longitude\ \(deg\) " + dat_fname
os.system(bashCmd + " >> " + tmp_fname)

# Need to use: subprocess library (e.g. subprocess.call())

# Extracting longitude and latitude data from RRI ephemeris.
import h5py
f = h5py.File(dat_fname)
geog_longs = f['CASSIOPE Ephemeris']['Geographic Longitude (deg)'].value
geog_lats  = f['CASSIOPE Ephemeris']['Geographic Latitude (deg)'].value

print "First Geographic Longitude: " + str(geog_longs[0])
print "First Geographic Latitude: " + str(geog_lats[0])

# Running Ashton's code for fun
import range_cells
#beams,gates = range_cells.findRangeCell([geog_lats[0]],[geog_longs[0]]) #for now just the first lat and long
start = timeit.default_timer()
beams,gates = range_cells.findRangeCell([0],[0],myFov=myFov)
end = timeit.default_timer()
print "Ashton code for one lat/lon pair took: " + str(end-start) + " seconds."

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

# Time the brute-force method
start = timeit.default_timer()

""" # Using script to go through 11 lat/lon pairs with Angeline Code
# Now, spin through each operational radar and test if it reaches (0,0)
nw = pydarn.radar.network()
results = dict()
lat_subset = geog_lats[0:10]
lon_subset = geog_longs[0:10]
for rad in nw.radars:
    if rad.status == 1:
        bm,gt,view=rgCoords.pos_to_rg(lat_subset,lon_subset,alt=300.0,
coords="geo",dtime=datetime(2016,04,01),frang=180.0,rsep=45.0,rad_id=rad.id) 
        results[rad.name]=(bm,gt,view)
"""

# Using script to go through 11 lat/lon pairs with Ashton code
nw = pydarn.radar.network()
results = dict()
lat_subset = geog_lats[0:10]
lon_subset = geog_longs[0:10]
for rad in nw.radars:
    if rad.status == 1:
        fov_f = davitpy.pydarn.radar.radFov.fov(site=rad.sites[0],altitude=300.,
            model='IS',coords='geo',ngates=75,fov_dir='front')
        bm_f,gt_f=range_cells.findRangeCell(lat_subset,lon_subset,fov_f)

        fov_b = davitpy.pydarn.radar.radFov.fov(site=rad.sites[0],altitude=300.,
            model='IS',coords='geo',ngates=75,fov_dir='back')
        bm_b,gt_b= range_cells.findRangeCell(lat_subset,lon_subset,fov_b)
        results[rad.name]=(bm_f,gt_f,bm_b,gt_b)

 
end = timeit.default_timer()
print "Time required to compute detailed intersections by brute force: " + \
    str(end - start) + " seconds."

# Takes 49 seconds doing the full brute-force calculations with input of (0,0)
# show results. For (0,0) input, no radar should reach it (maybe azores could).
for item in results.items():
    print item
    



"""
---------------------------------- PART 3 --------------------------------------
                      GRABBING RELEVANT SUPERDARN DATA             
--------------------------------------------------------------------------------
Having determined radars, beams, gates, grab the relevant data files (errlog &
perhaps more?) to include here.

Generate or grab a plot showing the satellite's track overlaid vs the superdarn
beams?


"""






