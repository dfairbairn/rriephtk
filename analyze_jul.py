"""
file: 'analyze_jul.py'
date: July 8 2016
author: David Fairbairn

Script with some analysis routines performed on the July 8th 2014 RRI data.
The goal of this analysis is to compare the received signal power from the 
various beams and infer which beam was operating at exactly 01:15:15 UTC for
the Saskatoon SuperDARN radar (we suspect the 1st or the 4th).

"""

import h5py
import numpy as np
import matplotlib.pyplot as plt

f = h5py.File("./data/RRI_20140708_011514_011741_lv1_v2.h5")

rri_data = f['RRI Data']
rri_settings=  f['RRI Settings']
radio_monopole1 = rri_data['Radio Data Monopole 1 (mV)']
radio_monopole2 = rri_data['Radio Data Monopole 2 (mV)']
radio_monopole3 = rri_data['Radio Data Monopole 3 (mV)']
radio_monopole4 = rri_data['Radio Data Monopole 4 (mV)']
ephem_times = f['CASSIOPE Ephemeris']['Ephemeris MET (seconds since May 24, 1968)'].value

# To get the ephems_to_datetime() method
from script_utils import *

rdm1 = radio_monopole1.value
rdm2 = radio_monopole2.value
rdm3 = radio_monopole3.value
rdm4 = radio_monopole4.value

""" 
# Wrapping up each of rdm.T[i] into one array as below seems to give wrong answers. 
# What the eff is the reason for these differing datasets within rdm1, rdm2, rdm3, rdm4?
allrdm1 = np.array(())
allrdm4 = allrdm3 = allrdm2 = allrdm1 # These are just copies, not the same pointer
for i in range(29):
    allrdm1 = np.concatenate([allrdm1,rdm1.T[i]])
    allrdm2 = np.concatenate([allrdm2,rdm2.T[i]])
    allrdm3 = np.concatenate([allrdm3,rdm3.T[i]])
    allrdm4 = np.concatenate([allrdm4,rdm4.T[i]])
"""
rdms = rdm4.T[1] # TBH I don't know what purpose is served by rdm1.T[i] for i = 1,2.. 
size = np.size(rdms)
time0 = ephem_times[0]

samprate = 2155.  #2083.34464    #2169.3356 #/62500.33933 

times_ephem = np.array(range(size))/samprate + time0 
times = ephems_to_datetime(times_ephem)

# For the purposes of plotting this stuff with each beam as a different color,
# I need to partition the dataset 13 different ways. A dictionary will store
# the 13 different subsets for both the time values and rdm values.
time_dict = dict()
rdm1_dict = dict()
for i in range(13):
    time_dict[i] = np.array([])
    rdm1_dict[i] = np.array([])

# Construct the time sets
# RRI doesn't start recording at the beginning of a beam transmission, so I 
# determined what the offset is by which beam demarcations are situated.
offset = 0 #int(0.2*samprate) #416
time_dict[12] = np.concatenate([time_dict[12],times_ephem[0:offset]])
relevant_rdms = rdms[0:offset]
rdm1_dict[12] = np.concatenate([rdm1_dict[12],relevant_rdms[0:offset]])
for i in range(146):
    lim1 = int(i*samprate + offset)
    lim2 = int((i+1)*samprate + offset) if i < 145 else 316722
    color_num = i % 13
    relevant_times = times_ephem[lim1:lim2]
    relevant_rdms = rdms[lim1:lim2]
    time_dict[color_num] = np.concatenate([time_dict[color_num],relevant_times])
    rdm1_dict[color_num] = np.concatenate([rdm1_dict[color_num],relevant_rdms])


#color_list = [(0.4 + x/5.0,0.4 + y/5.0, 0.4 + z/5.0) for x in range(3) for y in range(3) for z in range(3)]
# Specifically creating a color list with colors that contrast eachother
color_list = ['#808080','#FF0000','#FFFF00','#00FF00','#00FFFF','#0000FF',  
    '#FF00FF','#000000','#800000','#808000','#008000','#008080','#000080']

for i in range(13):
    plt.plot(time_dict[i],abs(rdm1_dict[i]),color=color_list[i])
plt.savefig("./data/output/20140708_sr" + str(int(samprate)) + "_ofs" + str(offset)+".jpg")
plt.show()
