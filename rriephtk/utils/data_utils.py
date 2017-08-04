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

import inspect
print("inspect.stack output: {0}".format(inspect.stack()[0][1]))
RRITK_ROOT = "../.."
RRITK_DATA = RRITK_ROOT + "/data"
RRITK_OUTPUT = RRITK_DATA + "/output"
RRITK_REMOTE = RRITK_DATA + "/remote"
DEFAULT_RRI = "RRI_20160418_222759_223156_lv1_v2.h5"

class FileLineWrapper(object):
    """
    A wrapper class for 'file' which tracks line numbers.

    **ATTRIBUTES**
        f (file object): the file object being wrapped.
        line (int): tracks current line number
        
    ** 0 indexed lines in the file! ** 
    ** Can use this just by going:
            flw = FileLineWrapper(open(fname,'r'))
    """
    def __init__(self, f):
        """ Constructor for filewrapper object. """
        self.f = f
        self.line = 1 
        
        # To allow skipping directly to lines
        self.line_offs = []
        offset = 0
        self.f.seek(0)
        for line in f:
            #print len(self.line_offs),line
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
        self.f.seek(self.line_offs[line_num-1])

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
    # TODO: Try-Excepts in here, use of subprocess more, 
    # Test if directory structure with ./data, ./data/output, ./data/remote exist
    if not os.path.isdir(RRITK_DATA):
        os.system('mkdir ' + RRITK_DATA)
        print "Creating script directory 'data/'."
    if not os.path.isdir(RRITK_OUTPUT):
        os.system('mkdir ' + RRITK_OUTPUT)
        print "Creating script directory 'data/output/'."
    if not os.path.isdir(RRITK_REMOTE):
        os.system('mkdir ' + RRITK_REMOTE)
        print "Creating script directory 'data/remote/'."

    # Check which data file is to be used based on command-line arguments
    import sys
    if sys.argv.__len__() == 2 and isinstance(sys.argv[1], str):
        dat_fname = sys.argv[1]
    else:
        print "No RRI file specified - going with default..."
        dat_fname = RRITK_DATA + "/" + DEFAULT_RRI 
    if not os.path.exists(dat_fname):
        print "No RRI file by that name. Exitting."
        return None, None

    # Check to see if we're on Maxwell.  If not, mount Maxwell SuperDARN data remotely
    if subprocess.check_output(["hostname"]) == "maxwell\n":
        print "Running on Maxwell: no remote mounting is necessary."
        data_path = "/data/" # If on Maxwell, don't do mounting
    else:
        # The sshfs tool is used to mount Maxwell's data directory locally.
        # 
        # First, the script unmounts anything currently already mounted in the mounting
        # directory (using os.system so that if nothing is there it doesnt crash the 
        # script).
        os.system("fusermount -uq " + RRITK_REMOTE)
        print "Accessing data files on maxwell, enter your password: "
        output = subprocess.check_output(["sshfs", 
                        "fairbairn@maxwell.usask.ca:/data/", RRITK_REMOTE])
        data_path = RRITK_REMOTE + "/" # Set the data path to where we just mounted Maxwell data
    return data_path, dat_fname

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
    # Wasted processing when other functions call this but don't receive the converted times array
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

def get_ottawa_data(date_string):
    """
    Function making it convenient in interactive mode or in scripts to select 
    data (including the index/seconds into the pass at which the Faraday 
    rotation reversal occurs).

    Source: Inspection of Plots of Orientation Angle

    """
    if isinstance(date_string, type(None)): date_string="20160418"

    # **** CHOOSE ONE OF THESE RRI FILES THEN RUN THE SCRIPT ****
    if "20160418"==date_string:
        fname = RRITK_DATA + "/RRI_20160418_222759_223156_lv1_v2.h5" #18th
        index_reversal = 167 #for 18th # Stay
    elif "20160419"==date_string:
        #index_reversal = 178 #for 19th # Could go -1
        index_reversal = 177
        fname = RRITK_DATA + "/RRI_20160419_220939_221336_lv1_v2.h5" #19th
    elif "20160420"==date_string:
        index_reversal = 213 #for 20th # Stay
        fname = RRITK_DATA + "/RRI_20160420_215117_215514_lv1_v2.h5" #20th
    elif "20160421"==date_string:
        #index_reversal = 205 #?? for 21st? # Could go up +5
        index_reversal = 210
        fname = RRITK_DATA + "/RRI_20160421_213255_213652_lv1_v2.h5" #21st
    elif "20160422"==date_string:
        #index_reversal = 222 #for 22nd # Could go down -1 
        index_reversal = 221
        fname = RRITK_DATA + "/RRI_20160422_211435_211832_lv1_v2.h5" #22nd
    else:
        print("Invalid input date.")
        return None    
    return fname,index_reversal

def get_ottawa_data2(date_string):
    """
    Does same as get_ottawa_data() only this time, it also includes the index
    of the reversal in ellipticity angle (in addition to the orientation angle
    flipping re: Faraday rotation).

    Param: date_string - formatted like "20160422"

    Returns: fname, rev_idx_orientation, rev_idx_ellipt

    Source: Inspection of Plots.
    (Glenn's inspection of ellipticity reversal, both of ours for orientation reversal)

    """
    if isinstance(date_string, type(None)): date_string="20160418"

    # **** CHOOSE ONE OF THESE RRI FILES THEN RUN THE SCRIPT ****
    if "20160418"==date_string:
        fname = RRITK_DATA + "/RRI_20160418_222759_223156_lv1_v2.h5" #18th
        orientation_rev = 168 #for 18th # Stay
        ellipt_rev = 167
    elif "20160419"==date_string:
        orientation_rev = 177
        ellipt_rev = 158
        fname = RRITK_DATA + "/RRI_20160419_220939_221336_lv1_v2.h5" #19th
    elif "20160420"==date_string:
        orientation_rev  = 213 #for 20th # Stay
        ellipt_rev = 133 # ??? it also does something near 213 though
        fname = RRITK_DATA + "/RRI_20160420_215117_215514_lv1_v2.h5" #20th
    elif "20160421"==date_string:
        orientation_rev  = 210
        ellipt_rev = 125
        fname = RRITK_DATA + "/RRI_20160421_213255_213652_lv1_v2.h5" #21st
    elif "20160422"==date_string:
        orientation_rev  = 221
        ellipt_rev = 110 
        fname = RRITK_DATA + "/RRI_20160422_211435_211832_lv1_v2.h5" #22nd
    else:
        print("Invalid input date.")
        return None    
    return fname,orientation_rev,ellipt_rev

def load_density_profile(fname):
    """
    Function for loading data from an electron densities data file.    

    ***PARAMS***
        fname [string]: name (including path) of file containing densities

    ***RETURNS***
        datarr [list of np.ndarray of floats] arrays corresponding to different
                 latitudes for e- density at various altitudes.
        datlats [list of floats]: latitudes corresponding to each entry in datarr list
 
    * NOTE: 

    """
    import pandas as pd
    import numpy as np
    dat = pd.read_csv(fname)
    datarr2 = ([(d[0].split()) for d in dat.as_matrix()])
    datarr = []
    datlats = []
    for i, drow in enumerate(datarr2):
        datarr2[i].pop(0)
        datlats.append(float(datarr2[i].pop(0)))
        datarr.append(np.array([ float(d) for d in drow ]))
    return datarr, datlats 

def load_rob_ephemeris(fname):
    """
    Use to grab the ephemeris points in one of Rob's simulation data files

    ***PARAMS***
        fname [string]: name (including path) of file containing ephemeris pts 

    ***RETURNS**
        Roblons [list of floats]: longitudes
        Roblats [list of floats]: latitudes 
        Robalts [list of floats]: altitudes

    * This function was written for parsing Rob Gillies' simulation's files. *

    """
    f = open(fname)
    Roblons = []
    Roblats = []
    Robalts = []
    for line in f:
        l = line.split()
        try:
            Roblons.append(float(l[1]))
            Roblats.append(float(l[2]))
            Robalts.append(float(l[3]))
        except:
            continue    
    return Roblons,Roblats,Robalts


def get_density(lon, lat, alt, datarr, datlats):
    """
    Function for taking array of electron densities and returning the right one
    for a particular location/altitude.

    (You give us the data and we'll do the work for you!)

    **These electron densities are in terms of m^3**


    ***PARAMS***
        lon [float]: longitude to check density at [deg]
        lat [float: latitude to check density at [deg]
        alt [float]: altitude to check density at [km]

    ***RETURNS***
        <density> [float]: electron number density from the densities  array [m^-3] 
    """
    import numpy as np
    if lat > np.max(datlats) or lat < np.min(datlats):
        #print("No data for that latitude - using closest latitude...")
        lat = np.max(datlats) if lat > np.max(datlats) else np.min(datlats)
    if alt > 559 or alt < 60:
        print("No data for that altitude")
        return -1
    for i,dlat in enumerate(datlats):
        alt_index = int(alt - 60)
        if lat < dlat:
            return datarr[i-1][alt_index]
        elif lat == dlat:        
            return datarr[i][alt_index]
    return 0

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
        return None
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


def list_mgf_files(path=(RRITK_DATA + '/mgf/')):
    """ 
    Returns a list of the mgf files that can simply be indexed so as to easily
    grab complicated filenames.
    """
    import os
    import re
    f_list = os.listdir(path)
    candidate_files = [ u for u in f_list if re.search("MGF_", u) ]
    return candidate_files

def exit_rri():
    """
    Easy shortcut in interactive mode to unmount the system and quit.

    ***PARAMS***
        - none -

    ***RETURNS***
        - none - 
    """
    import sys
    os.system("fusermount -uq " + RRITK_REMOTE)
    sys.exit()

# ----------------------------------------------------------------------------
#                            General Data Utilities
# ----------------------------------------------------------------------------

def ephems_to_datetime(ephem_times):
    """
    This function allows a whole array of ephemeris times to be conveniently
    converted in the same manner as with the ephem_to_datetime() function

    ** ARGS **
        ephem_times (numpy.array): array of Ephemeris MET (truncated JD times)

    ** RETURNS **
        times (list): list of datetime objects
    """
    import numpy as np
    if type(ephem_times)==list:
        ephem_times = np.array(ephem_times)
    if type(ephem_times)!=np.ndarray:
        print("Not an array")
        return None
 
    # i) Check the seconds between May 24 1968 (ephem MET) and Jan 1 1970 (neg number)
    t_off = subprocess.check_output(["date", "--date=1968-05-24 0:00:00", "+%s", "-u"])
    t_off = float(t_off.split("\n",1)[0]) # extract the integer value

    # ii) Do the math.
    ephtimes = ephem_times + t_off
    times = []
    for i in range(np.size(ephtimes)):
        times.append(dt.datetime.utcfromtimestamp(ephtimes[i]))
    return times

def ephem_to_datetime(ephem):
    """
    This function exists in order to conveniently convert the weird ephemeris
    time data format from the EPOP RRI instrument into a datetime object.

    *note*: makes use of the bash command 'date', so make sure you have it.

    ** ARGS **
        ephem (int): a large integer showing # of seconds since May 24 1968
                    (must be greater than 0, smaller than sys.maxint)

    ** RETURNS **
        dtime (datetime.datetime): a datetime object for the time given in ephem 
    """ 
    assert ephem >= 0
    import numbers
    assert isinstance(ephem, numbers.Number)

    # i) Check # seconds between May 24 1968 (ephem MET) and Jan 1 1970 (neg number)
    # -u parameter required to ensure offset is in UTC. +%s specifies output in seconds.
    t_off = subprocess.check_output(["date", "--date=1968-05-24 0:00:00", "+%s", "-u"])
    t_off = float(t_off.split("\n",1)[0]) # extract the integer value
  
    # ii) Do the math. 
    dtime = dt.datetime.utcfromtimestamp(t_off + float(ephem))    
    return dtime

def get_fov_by_name(name,frontback='front'):
    """
    This function shortcuts the FOV creation process by handling assumptions I
    typically make when creating FOV objects.

    ** ARGS **
        name (string): the name (as DaVitPy understands them) of a SuperDARN 
            radar whose FOV is to be returned.
        [frontback] (string): indicates whether to take front or back FOV 
    ** RETURNS **
        fov (pydarn.radar.radFov.fov): a (front) FOV object for the SuperDARN
            site with the given name.
    """
    from davitpy import pydarn
    nw = pydarn.radar.network()
    rad = nw.getRadarByName(name)
    assert (rad != False) #if name is illegal, getRadarByName returns False
    site = pydarn.radar.site(radId = rad.id)
    fov = pydarn.radar.radFov.fov(site=site,altitude=300.0,model='IS',
                                  coords='geo',ngates=75,fov_dir=frontback)
    return fov

def two_pad(in_time):
    """ 
    Takes in a number of 1 or 2 digits, returns a string of two digits. 
    """
    assert isinstance(in_time,int)
    assert (in_time < 100 and in_time >= 0)
    return "0" + str(in_time) if str(in_time).__len__() == 1 else str(in_time)

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    from math import radians, cos, sin, asin, sqrt
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r

def update_progress(progress):
    """
    # update_progress() : Displays or updates a console progress bar
    Accepts a float between 0 and 1. Any int will be converted to a float.
    A value under 0 represents a 'halt'.
    A value at 1 or bigger represents 100%

    Code by Brian Khuu
    """
    import time, sys
    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()


"""
TESTING
"""
if __name__ == "__main__":
    data_path,dat_fname = initialize_data()

    # Testing FileLineWrapper
    f = open('./script_utils.py','r')
    ln = f.readline()
    f.readline()
    offs = f.tell()
    line = f.readline()

    fw = FileLineWrapper(f)
    assert(fw.line_offs.__len__() > 0)
    assert(fw.line_offs[2] == offs)
    fw.seekline(3) 
    assert(fw.readline() == line)

    linen1 = fw.line
    fw.readline()
    assert(linen1 + 1 == fw.line)
    fw.seekline(0)
    fw.readline()

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
   

    # TESTS OF GENERAL DATA UTILITIES 
    # Testing ephem_to_datetime():
    tst_dt1 = dt.datetime(1968,5,24)
    assert(ephem_to_datetime(0) == tst_dt1)
    tst_dt2 = dt.datetime(2014,7,8,1,15,9)
    num2 = 1455498909
    assert(ephem_to_datetime(num2) == tst_dt2)
 
    # Testing ephems_to_datetime():
    import numpy as np
    assert(ephems_to_datetime(np.array([0,num2])) == [tst_dt1, tst_dt2])

    # Testing two_pad():
    assert(two_pad(3) == "03")
    assert(two_pad(13) == "13")

    # Testing get_fov_by_name():
    from davitpy import pydarn
    assert isinstance(get_fov_by_name("Saskatoon"),pydarn.radar.radFov.fov)

    exit_rri()
