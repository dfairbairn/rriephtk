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

import logging

# Loading data from July 8 2014
f = h5py.File("./data/RRI_20140708_011514_011741_lv1_v2.h5")
rri_data = f['RRI Data']
rri_settings=  f['RRI Settings']
radio_monopole1 = rri_data['Radio Data Monopole 1 (mV)']
radio_monopole2 = rri_data['Radio Data Monopole 2 (mV)']
radio_monopole3 = rri_data['Radio Data Monopole 3 (mV)']
radio_monopole4 = rri_data['Radio Data Monopole 4 (mV)']
ephem_times = f['CASSIOPE Ephemeris']['Ephemeris MET (seconds since May 24, 1968)'].value

rdm1 = radio_monopole1.value
rdm2 = radio_monopole2.value
rdm3 = radio_monopole3.value
rdm4 = radio_monopole4.value

sz_samps = np.shape(rdm1)[0]
sz_sets = np.shape(rdm1)[1]

# *** PARAMETER TO SET ***  
# Set this to RDM1, RDM2, RDM3, or RDM4 to select desired monopole data.
desired_dataset = rdm4

# For some reason the data is split up weirdly into 29 sets. In case I want 
# to look at all of it at once, I combine it together into allrdm below.
allrdm = np.zeros(sz_samps*sz_sets)
for i in range(sz_samps): 
    for j in range(sz_sets): # Sift together adjacent values in different sets.
        allrdm[sz_sets*i+j] = (desired_dataset.T[j])[i]

# *** PARAMETER TO SET ***
# 'rdms' will be whichever data set we choose to analyze.
# Look at a subset (eg. rdm1.T[0]) or allrdm (might be each rdm1.T[i] combined)
rdms = allrdm  
size = np.size(rdms)
print 'Data set size: ' + str(size)

# Samprate, approximately is # samples divided by the 146 seconds for the data
samprate = 62479. # ~2169 or ~62910. Instead, 2155. and 62479. are good
time0 = ephem_times[0]
times_ephem = np.array(range(size))/samprate + time0 
print 'Calculated Sample rate: ' + str(samprate)

# For the purposes of plotting this stuff with each beam as a different color,
# I need to partition the dataset 13 different ways. A dictionary will store
# the 13 different subsets for both the time values and rdm values.
time_dict = dict()
rdm_dict = dict()
for i in range(13):
    time_dict[i] = np.array([])
    rdm_dict[i] = np.array([])

# Construct the time sets for each different beam.
#
# RRI won't necessarily start recording at the beginning of a beam transmission 
# so an offset may be required to ensure beam partitions line up with the gap
# between the recorded beam sequences.
offset = 0 # Turns out, for this data, no offset is needed!
time_dict[12] = np.concatenate([time_dict[12],times_ephem[0:offset]])
relevant_rdms = rdms[0:offset]
rdm_dict[12] = np.concatenate([rdm_dict[12],relevant_rdms[0:offset]])
for i in range(146):
    lim1 = int(i*samprate + offset)
    lim2 = int((i+1)*samprate + offset) if i < 145 else (size - 1)
    color_num = i % 13
    relevant_times = times_ephem[lim1:lim2]
    relevant_rdms = rdms[lim1:lim2]
    time_dict[color_num] = np.concatenate([time_dict[color_num],relevant_times])
    rdm_dict[color_num] = np.concatenate([rdm_dict[color_num],relevant_rdms])

# Specifically creating a color list with colors that contrast eachother
color_list = ['#808080','#FF0000','#FFFF00','#00FF00','#00FFFF','#0000FF',  
    '#FF00FF','#000000','#800000','#808000','#008000','#008080','#000080']

# Finally, plotting each data sequence.
for i in range(13):
    plt.plot(time_dict[i],abs(rdm_dict[i]),color=color_list[i]) 
plt.savefig("./data/output/20140708_sr" + str(int(samprate)) + "_ofs" + str(offset)+".jpg")
plt.show()


