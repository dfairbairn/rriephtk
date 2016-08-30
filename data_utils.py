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

class FileLineWrapper(object):
    """
    A wrapper class for 'file' which tracks line numbers.

    **ATTRIBUTES**
        f (file object): the file object being wrapped.
        line (int): tracks current line number
        
    ** 0 indexed lines in the file! ** 
    """
    def __init__(self, f):
        """ Constructor for filewrapper object. """
        self.f = f
        self.line = 0
        
        # To allow skipping directly to lines
        self.line_offs = []
        offset = 0
        self.f.seek(0)
        for line in f:
            #print line
            self.line_offs.append(offset)
            offset += line.__len__()
        self.f.seek(0)
    def close(self):
        """ Closer function for filewrapper object. """
        return self.f.close()
    def readline(self):
        """ 
        This wrapper adds a lot of overhead with the if-statement, but has
        nicer functionality by checking if the file is empty before increasing line#...
        """
        ln = self.f.readline()
        if ln!="":
            self.line += 1
        return ln
    def seekline(self,line_num):
        """ Go to the nth line in the file. """ 
        self.line = line_num
        self.f.seek(self.line_offs[line_num])

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

def get_line_in_file(fl, srch_str):
    """
    A function taking a search string and a filelinewrapper obj as parameters, 
    returning the line number of the first occurrence of the search string in the file.

    *** PARAMS ***
    fl (FileLineWrapper):   a FileLineWrapper object instance to e searched
    srch_str (string):      the string of interest to find  

    *** RETURNS ***
    line_num (int): the line_num (of a FileLineWrapper object) where the search
                    string first occurs
    """
    assert isinstance(fl, FileLineWrapper)
    buffer_line = fl.line
    fl.line = 0
    found = False
    
    while 1:
        line_num = fl.line
        ln = fl.readline()
        if ln=="":
            print "End of file"
            break
        elif ln.find(srch_str) != -1:
            found = True
            break
        #if fl.line%1000 == 0:
        #    print fl.line
    fl.seekline(buffer_line)
    if False==found:
        return -1
    return line_num

def exit_rri():
    """
    Easy shortcut in interactive mode to unmount the system and quit.
    """
    import sys
    os.system("fusermount -uq ./data/remote/")
    sys.exit()

"""

TESTING

"""
if __name__ == "__main__":
    data_path,dat_fname = initialize_data()

    start = dt.datetime(2014,7,8,1,15,9)
    end = dt.datetime(2014,7,8,1,17,30)
    #start = dt.datetime(2014,7,8,1,0,0)
    #end = dt.datetime(2014,7,8,1,1,2,0)

    file_errl = open_errlog(data_path, 'sas', start)
    file_stamps = open_tstamps(data_path,start)

    # TODO: test get_line_in_file()??
    #print get_line_in_file(file_stamps, '01:15')
    #print get_line_in_file(file_stamps, '01:15')

    #TODO: automated tests on some of the new data access functions??? 
    exit_rri()
