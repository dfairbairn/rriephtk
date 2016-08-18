"""
file: 'data_utils.py'
author: David Fairbairn
date: July 2016

With a couple of different scripts that require some similar initial checks,
each of which filled with a lot of different work being done, it became clear
an initialization file would be required. This will help allow me to ensure
that portability is all good.

"""

import os
import subprocess

import datetime as dt
import logging 

from script_utils import * # For two_pad... easily replaceable functionality

def initialize_data():
    """
    Used by RRI scripts to ensure the expected directory structure and data
    exist and can be used.

    ** PARAMS ** 
        - none - 

    ** RETURNS **
        data_path (string):     A string containing the path to the root Maxwell
                                data, either remote or otherwise.
        data_fname (string):    A string containing the path and filename for the
                                RRI data e.g. ./data/RRI_20150402_032244_033241_lv1_v2.h5 
    """
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

    # Check which data file is to be used based on command-line arguments
    # TODO: allow initialize_data() to take an RRI file arg, resolve conflict of command-line and functional args.
    import sys
    if sys.argv.__len__() == 2 and isinstance(sys.argv[1], str):
        dat_fname = sys.argv[1]
    else:
        print "No RRI file specified - going with default..."
        dat_fname = "./data/RRI_20160418_222759_223156_lv1_v2.h5" # An RRI data file
    
    if not os.path.exists(dat_fname):
        print "No RRI file by that name. Exitting."
        exit()

    # Check to see if we're running on Maxwell, and if not, mount Maxwell SuperDARN data remotely
    if subprocess.check_output(["hostname"]) == "maxwell\n":
        print "Running on Maxwell: no remote mounting is necessary."
        data_path = "/data/" # If on Maxwell, don't do mounting
    else:
        # The sshfs tool is used to mount Maxwell's data directory locally.
        # 
        # First, the script unmounts anything currently already mounted in the mounting
        # directory (using os.system so that if nothing is there it doesnt crash the 
        # script).
        os.system("fusermount -uq ./data/remote/")
        print "Accessing data files on maxwell, enter your password: "
        output = subprocess.check_output(["sshfs", "fairbairn@maxwell.usask.ca:/data/","./data/remote"])
        data_path = "./data/remote/" # Set the data path to where we just mounted Maxwell data
 
    return data_path,dat_fname

def get_rri_ephemeris(dat_fname):
    """
    A function which returns the important ephemeris data, given an RRI h5 file.

    *** PARAMS *** 
    dat_fname (string): string giving the relative path and filename for the RRI h5 file

    *** RETURNS ***
    glons (float): Geographic Longitude (degrees)
    glats (float): Geographic Latitude (degrees) 
    alts (float): Altitude (mk)
    etimes (float): Ephemeris MET/Truncated JD time (seconds since May 24 1968)
    """
    import h5py
    f = h5py.File(dat_fname)
    geog_longs = f['CASSIOPE Ephemeris']['Geographic Longitude (deg)'].value
    geog_lats  = f['CASSIOPE Ephemeris']['Geographic Latitude (deg)'].value
    ephem_times = f['CASSIOPE Ephemeris']['Ephemeris MET (seconds since May 24, 1968)'].value
    alts = f['CASSIOPE Ephemeris']['Altitude (km)'].value
    return geog_longs,geog_lats,alts,ephem_times

def get_rri_ephemeris_full(dat_fname):
    """
    A similar accessor function for RRI ephemeris from the relative path+name 
    of an RRI file, but this time also grabbing the MLT, MLAT, MLONG data.

    *** PARAMS ***
    dat_fname (string): string giving the relative path and filename for the RRI h5 file

    *** RETURNS ***
    glons (float): Geographic Longitude (degrees)
    glats (float): Geographic Latitude (degrees) 
    alts (float): Altitude (mk)
    etimes (float): Ephemeris MET/Truncated JD time (seconds since May 24 1968)
    mlon (float): Magnetic longitude (degrees)
    mlat (float): Magnetic latitude (degrees)
    mlts (float): Magnetic Local Time (hr)
    pitch (float): pitch of CASSIOPE (deg)
    yaw (float): yaw of CASSIOPE (deg)
    roll (float): roll of CASSIOPE (deg)
    """
    import h5py
    from davitpy.models import aacgm
    f = h5py.File(dat_fname)
    geog_longs = f['CASSIOPE Ephemeris']['Geographic Longitude (deg)'].value
    geog_lats  = f['CASSIOPE Ephemeris']['Geographic Latitude (deg)'].value
    ephem_times = f['CASSIOPE Ephemeris']['Ephemeris MET (seconds since May 24, 1968)'].value
    #TODO: Wasted processing when other functions call this but don't receive the converted times array
    times = ephems_to_datetime(ephem_times) 
    alts = f['CASSIOPE Ephemeris']['Altitude (km)'].value
    mlat = f['CASSIOPE Ephemeris']['Magnetic Latitude (deg)'].value
    mlon = f['CASSIOPE Ephemeris']['Magnetic Longitude (deg)'].value
    mlts = []
    for i in range(mlon.__len__()):
        dt = times[i]
        lone_mlon = mlon[i]
        mlts.append(aacgm.mltFromYmdhms(dt.year,dt.month,dt.day,dt.hour,dt.minute,dt.second,lone_mlon))
    pitch = f['CASSIOPE Ephemeris']['Pitch (deg)'].value
    yaw = f['CASSIOPE Ephemeris']['Yaw (deg)'].value
    roll = f['CASSIOPE Ephemeris']['Roll (deg)'].value
    return geog_longs,geog_lats,alts,ephem_times,mlon,mlat,mlts,pitch,yaw,roll


def get_hdf5(dat_fname):
    """
    A function which mostly just wraps the use of the h5py library.
    """
    import h5py
    f = h5py.File(dat_fname)
    return f

def open_errlog(data_path, rcode, date):
    """
    A function that can make opening errlog files a re-usable action?

    ** PARAMS **
        rcode (string):
        date (datetime object): the date of interest
        data_path (string): string containing the (possibly relative) path to the errlog data

    ** RETURNS **
        file_errlog (filelinewrapper object??):
    """
    import bz2
    fpath = data_path + rcode + "_errlog/" #/path/to/data/ + sas + _errlog/
    fname_reg = str(date.year) + str(two_pad(date.month)) + str(two_pad(date.day)) + "." + rcode + ".errlog"
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
    return FileLineWrapper(file_errl)

def open_tstamps(data_path, date):
    """
    A function that makes opening timestampdata for Saskatoon a reusable action.

    ** PARAMS **
        date (datetime object): the date of interest
        data_path (string): contains the (possibly relative) path to the greater data folder

    ** RETURNS **
        file_timestampdata (filelinewrapper object??): 
    """
    import bz2
    fname = str(date.year) + str(two_pad(date.month)) + str(two_pad(date.day)) \
     + "." + str(two_pad(date.hour)) + "00" + ".timestampdata.bz2"
    file_stamps = bz2.BZ2File(data_path + "epop/" + fname)
    return FileLineWrapper(file_stamps)

def get_stamp_pulses(file_stamps):
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
    return pulses

def get_diffs(pulses):
    """ Returns a list of time differences between the pulse times given, corresponding 1 - for - 1 """
    diffs = []
    for i in range(pulses.__len__() - 1):
        p = pulses[i]
        p_nxt = pulses[i+1]
        diff = p_nxt-p
        diffs.append(p_nxt-p)
    return diffs

def determine_offset(pulse_times_a, pulse_seqs_a, pulse_times_b, pulse_seqs_b):
    """
    This function determines a discrete offset by which the first list of pulse
    sequences can be shifted so as to make the pulses with the same index in 
    each list be most similar.
    e.g.
    [3,4,3,2]

    """
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
    print overlap_scores
    print overlap_shifts
    return best_overlap_index, 0.5  
    #TODO: Confidence/quality of answer???

def evaluate_difference(lst_a, lst_b):
    """
    Determines how much overlap there is between the two input lists.
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

"""

TESTING

"""
if __name__ == "__main__":
    data_path,dat_fname = initialize_data()
    file_errl = open_errlog(data_path, 'sas', dt.datetime(2014,7,8))
    file_tstamps = open_tstamps(data_path, dt.datetime(2014,7,8))
    lista = [7,8,8,8,8,7,8,8,8,8,7,8,8,8,8,7]
    listb = [8,8,8,7,8,8,8,8,7,8,8,8,8,7] 
    print "List A: " + str(lista)
    print "List B: " + str(listb)
    score = evaluate_difference(lista,listb)
    print "Difference evaluation: " + str(score)
    shift = determine_shift_offset(lista,listb)
   
    #TODO: automated tests on some of the new data access functions??? 
