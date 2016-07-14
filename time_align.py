"""
file: 'time_align.py'
author: David Fairbairn
date: 27th June 2016

The need for a script that looks at timestampdata (currently only relevant for
the Saskatoon SuperDARN radar) to the errlog files' erroneous timestamps
compelled me to write this script.

This script approaches the problem by identifying the times during which 7 and 
8 pulse sequences occur, and using a pattern of these (e.g. 7,8,7,7,7,8,..) to
find when the pattern begins in the ERRLOG file and when it begins in the 
reliable timestampdata file. The mapping of corresponding pulse sequences in 
each file allows us to deduce correct times for the errlog data. 

Likely due to running old software (the QNX operating system) running on new 
hardware at the radar site, the main SuperDARN system at a few locations 
undergoes frequent and unpredictable timing corrections (e.g. every 5 minutes 
on average, discrete corrections that average about 0.5 seconds. 


"""


import datetime as dt
import numpy as np

import os
import subprocess
import logging

from script_utils import *

import bz2

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

# The sshfs tool is used to mount Maxwell's data directory locally.
# 
# First, the script unmounts anything currently already mounted in the mounting
# directory (using os.system so that if nothing is there it doesnt crash the 
# script).
os.system("fusermount -uq ./data/remote/")
print "Accessing data files on maxwell, enter your password: "
output = subprocess.check_output(["sshfs", "fairbairn@maxwell.usask.ca:/data/","./data/remote"])

# TODO:Which time should we look at? 
# Currently: we will default to looking at 2014-07-08 0700h
date = dt.datetime(2014,7,8,1) # THIS DAY ONLY HAS 8 PULSE SEQUENCES ? ALSO: 0100 or 0700???
#date = dt.datetime(2015,4,2,3) # THIS DAY SHOULD HAVE BOTH. ONLY HAS BOTH AT 0300h (not 0900h)
#date = dt.datetime(2014,8,7,1)

# Open the Timestamp data
fname = str(date.year)+  str(two_pad(date.month)) + str(two_pad(date.day)) \
+ "." + str(two_pad(date.hour))+ str(two_pad(date.minute)) + ".timestampdata.bz2"

file_stamps = bz2.BZ2File("./data/remote/epop/" + fname)


# Open the Saskatoon Errlog
rcode = 'sas' 
fpath = "./data/remote/" + rcode + "_errlog/"
fname_reg = str(date.year) + str(two_pad(date.month)) + str(two_pad(date.day)) \
                                                 + "." + rcode + ".errlog"
fname_bz2 = fname_reg + ".bz2"

if os.path.exists(fpath + fname_bz2):
    file_errl = bz2.BZ2File(fpath + fname_bz2)
    logging.debug("Uses bz2 errlog file.")
elif os.path.exists(fpath + fname_reg):
    logging.debug("Uses non-bz2 errlog file.")
    file_errl = open(fpath + fname_reg)
else:
    logging.error("No ERRLOG file found!") 
    exit()

# Reading Timestamp data, acquiring timing differences
end = False
pulses = []
while end != True:
    ln = file_stamps.readline()
    if ln == '':
        end = True
    elif ln.find("SEC") != -1:
        time = float((ln.split(" = ")[1]).split("\n")[0])
        pulses.append(time)

diffs = []
for i in range(pulses.__len__() - 1):
    p = pulses[i]
    p_nxt = pulses[i+1]
    diff = p_nxt-p
    diffs.append(p_nxt-p)

summary = []
pulse_seq = []
i = 0
minutes_count = 0
while i < diffs.__len__():
    d1 = diffs[i] 

    # Implemented a hack to notice minute-to-minute transitions, note them in summary file
    if d1 < 0: # This may screw up 7 or 8 pulse identification across transitions 
        d1 = d1 + 60.0 
        minutes_count = minutes_count + 1
        summary.append("\n" + two_pad(date.hour) + ":" + two_pad(minutes_count) + "\n")
    
    if i < diffs.__len__() - 6:
        d2 = diffs[i+1]
        d3 = diffs[i+2]
        d4 = diffs[i+3]
        d5 = diffs[i+4]
        d6 = diffs[i+5]
        d7 = diffs[i+6]

        c1 = np.around(d1,decimals=4) == 0.0210
        c2 = np.around(d2,decimals=4) == 0.012
        c3 = np.around(d3,decimals=4) == 0.003
        c4 = np.around(d4,decimals=4) == 0.0045
        c5 = np.around(d5,decimals=4) == 0.006
        c6 = np.around(d6,decimals=4) == 0.0165
        c7 = np.around(d7,decimals=4) == 0.0015
        
        b1 = np.around(d1,decimals=4) == 0.0216
        b2 = np.around(d2,decimals=4) == 0.0072
        b3 = np.around(d3,decimals=4) == 0.0192
        b4 = np.around(d4,decimals=4) == 0.0048
        b5 = np.around(d5,decimals=4) == 0.0096
        b6 = np.around(d6,decimals=4) == 0.0024

        if c1 and c2 and c3 and c4 and c5 and c6 and c7:
            #print "8 Pulse Sequence"
            summary.append("8 Pulse Sequence")
            pulse_seq.append("8")
            i = i+7
        elif b1 and b2 and b3 and b4 and b5 and b6:
            #print "7 pulse sequence"
            summary.append("7 Pulse Sequence")
            pulse_seq.append("7")
            i = i+6
        else:
            summary.append(str(d1))
            i = i + 1
    else:
        #print d
        summary.append(str(d1))
        i = i  + 1


f = open("./data/output/timestamp_sum_"+str(date.year)+str(two_pad(date.month)) + str(two_pad(date.day)) + str(two_pad(date.hour)) + ".dat",'w+')
f2 = open("./data/output/timestamp_pulses_"+str(date.year)+str(two_pad(date.month)) + str(two_pad(date.day)) + str(two_pad(date.hour)) + ".dat",'w+')

for s in summary:
    f.write(s + "\n")

for p in pulse_seq:
    f2.write(p + "\n")

# Reading the ERRLOG data!

#file_errl.readline()
ln = file_errl.readline()

# How much of errlog to analyze for this?

# the main logic:
"""
if ln.find("Number of Sequences [") != -1:
    lp = (ln.split(" : ")[2]).split("[")[1]
    pulse = int(lp.split("]")[0])
    numof = int((lp.split(" ")[1]).split("\n")[0])
""" 

# Unmount the Maxwell remote mount.
os.system("fusermount -uq ./data/remote/")
