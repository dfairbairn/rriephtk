"""
File: 'script-utils.py'
Description: Contains utilities used in the RRI Conjunction Finder script.
    These include functions to more conveniently grab FOV data from Davitpy,
    plotting functions that make use of Matplotlib, etc.
Author: David Fairbairn
Date: June 2016

"""
import matplotlib.pyplot as plt
import subprocess 
import numpy as np 
import datetime as dt #A wrapper class for 'file' which tracks line numbers. 

from davitpy.utils import plotUtils
def ephems_to_datetime(ephem_times):
    """
    This function allows a whole array of ephemeris times to be conveniently
    converted in the same manner as with the ephem_to_datetime() function

    ** ARGS **
        ephem_times (numpy.array): array of Ephemeris MET (truncated JD times)

    ** RETURNS **
        times (list): list of datetime objects
    """
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

    # Testing ephem_to_datetime():
    tst_dt1 = dt.datetime(1968,5,24)
    assert(ephem_to_datetime(0) == tst_dt1)
    tst_dt2 = dt.datetime(2014,7,8,1,15,9)
    num2 = 1455498909
    assert(ephem_to_datetime(num2) == tst_dt2)
 
    # Testing ephems_to_datetime():
    assert(ephems_to_datetime(np.array([0,num2])) == [tst_dt1, tst_dt2])

    # Testing two_pad():
    assert(two_pad(3) == "03")
    assert(two_pad(13) == "13")

    # Testing get_fov_by_name():
    from davitpy import pydarn
    assert isinstance(get_fov_by_name("Saskatoon"),pydarn.radar.radFov.fov)
