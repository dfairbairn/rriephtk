"""
File: 'beams_gates_usage.py'
Code from checking out Ashton and Angeline's beam/gate determiners
Author: David Fairbairn
Date: May 2016
"""

import davitpy
from davitpy import pydarn
from datetime import datetime
import numpy as np

dat_fname = "/home/david/pyth_stuff/script/data/RRI_20160401_072714_073111_lv1_v1.h5" # An RRI data file

#radars = davitpy.pydarn.network()

# Creating an FOV
site = pydarn.radar.site(code='rkn')
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
#beams,gates = range_cells.findRangeCell(lats,lons, myFov=fov)



# ANGELINE BURRELL COORDS->RG CODE
import rgCoords

dtime=datetime(2011,06,01)
# rg_gate: range gates, bmnum: beam numbers, fovflg: flag showing forward FOV or back lobe
# NOTE: rlats and rlons currently are actually lat/lon coords meant for 
#       checking FOV containment, not for being contained themselves. JUST A TEST
#rg_gate,bmnum,fovflg = rgCoords.pos_to_rg(rlats[7],rlons[7],coords='geo',alt=300.0,dtime=dtime,frang=180.0,rsep=45.0,rad_code='sas')





"""
Additional homeless code that doesn't belong in the script anymore but might be
useful for reference.
"""

from davitpy import pydarn
import numpy as np
import os

# A now-unnecessary bit of code to generate new data files (and not overwrite)
newfile = "data/output/tmp_ephem.dat" 
tmp_fname = newfile
i = 0
while (os.path.isfile(tmp_fname)):
    i+=1
    tmp_fname = newfile + str(i)


# This bash command is now superfluous, because we use h5py to get data 
bashCmd = "h5dump -d CASSIOPE\ Ephemeris/Geographic\ Longitude\ \(deg\) " + dat_fname
os.system(bashCmd + " >> " + tmp_fname)



import subprocess
# Ephemeris MET timing data is in the awkward format of 'seconds since May 24 1968'
# The current time can be converted to a similar number by doing as follows:
# i) Use the date command in Bashs option to give # of seconds since Epoch (1970 Jan 1st)
t1 = subprocess.check_output(["date", "+%s"])
t1 = int(t1.split("\n",1)[0]) # extract the integer value

# ii) Check the seconds between May 24 1968 and Jan 1 1970 (neg number)
t2 = subprocess.check_output(["date", "--date=1968-05-24 0:00:00", "+%s"])
t2 = int(t2.split("\n",1)[0]) # extract the integer value

# iii) Take the difference
t = t1-t2



