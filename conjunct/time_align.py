"""
file: 'time_align.py'
author: David Fairbairn
date: June 2016

The need for a script that looks at timestampdata (currently only relevant for
the Saskatoon SuperDARN radar) to the errlog files' erroneous timestamps
compelled me to write this script.

This script approaches the problem by identifying the times during which 7 and 
8 pulse sequences occur, and using a pattern of these (e.g. 7,8,8,8,7,8,..) to
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
import matplotlib.pyplot as plt

import os
import subprocess
import logging

from script_utils import *
from data_utils import * # for initialize_data(), open_tstamps(), open_errlog()

import bz2

# ======================= FUNCTIONS FOR TIME ALIGNMENT ========================

def get_stamp_pulses(file_stamps,start_time,end_time):
    """
    A function which 

    *** PARAMS ***
    file_stamps (file object): timestamps file e.g. maxwell:/data/epop/20160418.0100.timestampdata.bz2
    start_time (datetime obj): the start of the time period of interest for gathering pulse data
    end_time (datetime obj): the end of the time period of interest 

    *** RETURNS *** 
    pulse_times (list of strings): the total timestamps (as a string) for each pulse
    pulses (list of floats): the times (in seconds) at which a pulse goes out    
    
    """
    #TODO: GET CLOSER TO THE ACTUAL START < 10 s (rather than as with these params, up to 59 seconds away)
    strt_str = two_pad(start_time.hour) + ":" + two_pad(start_time.minute)
    #TODO: Find way to grab minute after the one of interest
    end_str = two_pad(end_time.hour) + ":" + two_pad(end_time.minute)

    startln = get_line_in_file(file_stamps,strt_str)
    endln = get_line_in_file(file_stamps,end_str)

    print "Start line for search string of " + strt_str + ": " + str(startln)
    print "End line for search string of " + end_str + ": " + str(endln)

    # Reading Timestamp data, acquiring timing differences

    end = False
    pulse_times = []
    pulses = []
    # Initialized hour/minute timestamp for edge cases where a pulse is read 
    # in before it has a corresponding hr/min
    hrtime = "--:--"
    file_stamps.seekline(startln)
    while end != True:
        ln = file_stamps.readline()
        if ln == '' or file_stamps.line > endln:
            print "End of file or reached end of search range."
            end = True
        elif ln.find("TIME") != -1:
            hrtime = (ln.split(" = ")[1]).split(" ")[0]
        elif ln.find("SEC") != -1:
            sectime = float((ln.split(" = ")[1]).split("\n")[0])
            if sectime < 10.0:
                time = hrtime + ":0" + str(round(sectime,5))
            else:
                time = hrtime + ":" + str(round(sectime,5))
            pulse_times.append(time)
            pulses.append(sectime)
    return (pulse_times, pulses)

def get_errl_pulses(f_errl, start, end):
    """
    Function to grab the pulses from the errlog file for the desired time
    interval as well as their general timestamps.


    *** PARAMS ***
    file_errl (FileLineWrapper obj): errl file e.g. maxwell:/data/sas_errlog/...
    start (datetime obj): the start of the time period of interest for gathering pulse data
    end (datetime obj): the end of the time period of interest 

    *** RETURNS ***


    """
    start_str = two_pad(start.hour) + ":" + two_pad(start.minute) + ":"
    end_str = two_pad(end.hour) + ":" + two_pad(end.minute) + ":"
    #TODO: Find way to grab minute after the one of interest
    #end_str = two_pad(end.hour) + ":" + two_pad(end.minute + 1) + ":"

    ln_start = get_line_in_file(f_errl, start_str)
    ln_end = get_line_in_file(f_errl, end_str)

    print "Start line for search string of " + start_str + ": " + str(ln_start)
    print "End line for search string of " + end_str + ": " + str(ln_end)

    end = False
    pulse_times = []
    pulses = []

    f_errl.seekline(ln_start)
    while end != True:
        ln = f_errl.readline()
        if ln.find("Number of sequences") != -1:
            #print "Found pulse sequence!"
            pulse,numof = parse_pulses(ln)
            ptime = parse_ptimes(ln)
            for i in range(numof):
                pulses.append(pulse)
                pulse_times.append(ptime)
        elif ln == '' or f_errl.line > ln_end:
            print "End of file or reached end of search range."
            end = True
    return (pulse_times, pulses)

def get_diffs(pulse_times, pulses):
    """ 
    Returns a list of time differences between the pulse times given, 
    corresponding 1 - for - 1 

    *** PARAMS ***
    pulse_time (list of strings): list of strings containing the exact time of the pulse.
    pulses (list of floats): list of floats of the exact time (sec) of the pulses.

    *** RETURNS ***
    diff_val ([] floats): time intervals between temporally adjacent pulses.
    pulse_times ([] strings): timestamp of beginning of each time interval
    """
    
    diffs = []
    for i in range(pulses.__len__() - 1):
        p = pulses[i]
        p_nxt = pulses[i+1]
        diff = p_nxt-p
        diffs.append(p_nxt-p)
    return pulse_times,diffs

def identify_sequences(pulse_times,diffs):
    """
    This function takes a list of time intervals between pulses (whose interval
    begins at the corresponding time in pulse_times, and picks out which series
    of intervals corresponds to a 7 or 8 pulse sequence.

    *** PARAMS ***
    diff_val ([] floats): time intervals between temporally adjacent pulses.
    pulse_times ([] strings): timestamp of beginning of each time interval

    *** RETURNS ***
    total ([] strings): list of every feature in diffs (possibly deprecated)
    sequence_times ([] strings): timestamp of beginning of each time interval
    sequences ([] strings): 7 vs 8 for which sequence occurred

    """
    total = []
    sequence_times = []
    sequences = []
    i = 0
    # minutes_count separately tracks the number of minute-to-minute transitions
    # that the sequence-identifier finds (which should hopefully match what's in pulse_times)
    minutes_count = 0
    while i < diffs.__len__():
        d1 = diffs[i] 
        t1 = pulse_times[i] 
        # Implemented a hack to notice minute-to-minute transitions, note them in summary file
        if d1 < 0: # This may screw up 7 or 8 pulse identification across transitions 
            d1 = d1 + 60.0 
            minutes_count = minutes_count + 1
            total.append("Minute Transition: " + str(minutes_count))
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
                total.append("8 Pulse Sequence")
                sequence_times.append(t1)
                sequences.append("8")
                i = i+7
            elif b1 and b2 and b3 and b4 and b5 and b6:
                #print "7 pulse sequence"
                total.append("7 Pulse Sequence")
                sequence_times.append(t1)
                sequences.append("7")
                i = i+6
            else:
                total.append(str(d1))
                i = i + 1
        else:
            #print d
            total.append(str(d1))
            i = i  + 1
    return total,sequence_times,sequences

def parse_pulses(ln):
    """ 
    A little mini function for taking a line in an errlog file and grabbing the pulse number.
    """
    if ln.find("Number of sequences") == -1:
        print "Line of text doesn't contain pulse information!"

    rem = ln.split("Number of sequences")[1]
    if rem.find("[") != -1:
        # Then this file *does* specify pulse sequences as it should
        pseq = (rem.split("[")[1]).split("]")[0]
        numof = int((rem.split(": ")[1]).split("\n")[0])
    else:
        # Then this is an older errlog where the pseqs are all 8 pulse sequences
        pseq = str(8)
        numof = int((rem.split(": ")[1]).split("\n")[0])
    return pseq,numof 

def parse_ptimes(ln):
    """
    Mini function for getting the time (in a string) from a line of an errlog file.
    """
    timestring = ln.split(" ")[3]
    return timestring

def determine_offset(pulse_times_a, pulse_seqs_a, pulse_times_b, pulse_seqs_b):
    """
    This function determines a discrete offset by which the first list of pulse
    sequences can be shifted so as to make the pulses with the same index in 
    each list be most similar.
    e.g.
    [3,4,3,2]

    """
    #TODO: Finish this function
    return -1, 0

def determine_shift_offset(lst_a, lst_b):
    """
    Determines the optimal shifting of one sequence with respect to the other
    so that the most list entries with the same indices are equal between the
    two lists.
    
    ** PARAMS **
        lst_a (list): the first list, which the function hopes is the smaller.
        lst_b (list): the second list.

    ** RETURNS **
        offset (integer): the offset from lst_a with respect to lst_b yielding
                        optimal matching between the two lists.
        quality (integer): 0 for poor confidence, 1 for strong confidence
    """
    assert isinstance(lst_a, list)
    assert isinstance(lst_b, list)
    lst_a_len =  lst_a.__len__()
    # Ensure lst_a is the shorter list
    if lst_a_len > lst_b.__len__():
        return determine_shift_offset(lst_b, lst_a) 
    
    # Loop through a reasonable different number of integer index shifts to try 
    # The first one we try should be no shift whatsoever, and if there's 100% 
    # overlap, don't bother trying anything else (???).
    #for i in range(lst_a_len):
    best_overlap = evaluate_difference(lst_a,lst_b)
    best_overlap_index = 0
    if best_overlap == 1.0:
        # Optimal overlap already
        confidence = 1
        return best_overlap_index, confidence
    #TODO: Should we collect several  ?

    overlap_scores = []
    overlap_shifts = []
    lst_b_len = lst_b.__len__()
    lst_b_fwd = lst_b_bck = lst_b
    for i in range(lst_a_len):
        lst_b_fwd = lst_b_fwd[1:lst_b_len] + [lst_b_fwd[0]] # rotate forward   
        overlap = evaluate_difference(lst_a, lst_b_fwd)
        if overlap > best_overlap:
            best_overlap = overlap
            # access indices start at 0, so subtract 1 to describe extent of shift
            best_overlap_index = -i - 1
        overlap_scores.append(overlap)
        overlap_shifts.append(-i-1)

        lst_b_bck = [lst_b_bck[-1]] + lst_b_bck[0:-1]        
        overlap = evaluate_difference(lst_a, lst_b_bck)
        if overlap > best_overlap:
            best_overlap = overlap
            # access indices start at 0, so subtract 1 to describe extent of shift
            best_overlap_index = i + 1
        overlap_scores.append(overlap)
        overlap_shifts.append(i+1)
    #print overlap_scores # Less output
    #print overlap_shifts
    return best_overlap_index
    #TODO: Confidence/quality of answer???

def evaluate_difference(lst_a, lst_b):
    """
    Determines how much overlap there is between the two input lists.

    Essentially a value function to maximize.
    """
    assert isinstance(lst_a, list)
    assert isinstance(lst_b, list)
    lst_a_len =  lst_a.__len__()
    if lst_a_len > lst_b.__len__():
        return evaluate_difference(lst_b, lst_a) 
    diff = 0.0
    for i in range(lst_a_len):
        if lst_a[i] != lst_b[i]:
            diff = diff + 1.0
    # The score will be an overlap percentage
    return (lst_a_len - diff)/lst_a_len

def visualize_list_difference(lst_a, lst_b):
    """

    """
    import collections
    l = lst_a
    y = []
    for i in range(lst_a.__len__()):
         d = collections.deque(l)
         d.rotate(-1)
         l = (np.array(d)).tolist()
         y.append(evaluate_difference(l,lst_b))
    return y


# ========================= TIME ALIGNMENT SCRIPT =============================

data_path,dat_fname = initialize_data() 

#start = dt.datetime(2014,7,8,1,15,9)
#end = dt.datetime(2014,7,8,1,17,30)
start = dt.datetime(2016,4,18,0,30,0)
end = dt.datetime(2016,4,18,0,33,0)

# Open the Timestamp data
file_stamps = open_tstamps(data_path, start)

# Open the Saskatoon Errlog
rcode = 'sas' # If we had another Timestamper, this could be an input parameter
file_errl = open_errlog(data_path, rcode, start)

# Reading Timestamp data, acquiring timing differences
stamp_ptimes,stamp_pulses = get_stamp_pulses(file_stamps, start, end)
stamp_dtimes,stamp_diffs = get_diffs(stamp_ptimes,stamp_pulses)
stamp_allpulses,stamp_seqtimes,stamp_pseqs = identify_sequences(stamp_dtimes,stamp_diffs) 

# Reading the ERRLOG data!
errl_seqtimes,errl_pseqs = get_errl_pulses(file_errl, start, end)    

print("\nNow defining custom lists lista and listb...")

lista = [7,8,8,8,8,7,8,8,8,8,7,8,8,8,8,7]
listb = [8,8,8,7,8,8,8,8,7,8,8,8,8,7] 
score = evaluate_difference(lista,listb)
print("'evaluate_difference' result on lista vs listb initially: {0}".format(score))
shift = determine_shift_offset(lista,listb)
print("determined shift offset: {0}".format(shift))

errl_pseqs_ab = errl_pseqs[:len(stamp_pseqs)]
errl_seqtimes_ab = errl_seqtimes[:len(stamp_seqtimes)]
print("\nNow printing timestamper and errlog sequences and times" + 
    " near the start and end to show their alignment...".format(len(errl_seqtimes)))
for i in np.arange(28,42):
    stamp_str = str(stamp_seqtimes[i]) + "\t" + str(stamp_pseqs[i])
    errl_str = "\t" + str(errl_pseqs_ab[i]) + "\t" + str(errl_seqtimes[i])
    print(stamp_str + errl_str)
print("\n\n")
for i in np.arange(60, 45, -1):
    stamp_str = str(stamp_seqtimes[-1-i]) + "\t" + str(stamp_pseqs[-1-i])
    errl_str = "\t" + str(errl_pseqs_ab[-1-i]) + "\t" + str(errl_seqtimes[-1-i])
    print(stamp_str + errl_str)

# TODO: figure out what I was going to do with these two lines:
indx_del = 70
chunks = int(errl_pseqs.__len__()/70.)

# Run the stats on the equal-length versions of this data
score = evaluate_difference(stamp_pseqs, errl_pseqs_ab)
print("'evaluate_difference' result on lista vs listb initially: {0}".format(score))
shift = determine_shift_offset(stamp_pseqs, errl_pseqs_ab)
print("determined shift offset: {0}".format(shift))


print("\nNow supposedly going to perform visualization of the list differences in terms of offset similarities")
y = visualize_list_difference(errl_pseqs_ab, stamp_pseqs)
plt.plot(100.0*np.array(y))
plt.xlabel('Discrete rotations of list 1 vs list 2')
plt.ylabel('Agreement (%)')
plt.show()
# Unmount the Maxwell remote mount.
#os.system("fusermount -uq ./data/remote/")

exit_rri()
