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
def ephems_to_datetime(ephem_times):
    """
    This function allows a whole array of ephemeris times to be conveniently
    converted in the same manner as with the ephem_to_datetime() function

    ** ARGS **
        ephem_times (numpy.array): array of Ephemeris MET (truncated JD times)

    ** RETURNS **
        times (list): list of datetime objects
    """
    assert isinstance(ephem_times, np.ndarray), "Not an array"
 
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

def plot_sat_ephemeris(date_string=None,lons=None,lats=None,ephtimes=None):
    """


    """
    if lons==None or lats==None or ephtimes==None:
        if date_string==None or not isinstance(date_string, type("e.g.")):
            print "Need to provide at least one parameter to go on!"
            return -1
        from ottawa_plots import *
        # TODO: Fix this!~!!!! We should have a general RRI file access method outside of 'get ottawa data'!!
        rri_fname,indx_rev = get_ottawa_data(date_string) 
        lons,lats,alts,ephtimes = get_rri_ephemeris(rri_fname)
        times=  ephems_to_datetime(ephtimes)
    
    # A different font for the legend etc. might be nice
    #fig = plt.figure()
    font = {'fontname':'Computer Modern'}
    m = plotUtils.mapObj(lat_0=np.mean(lats), lon_0=np.mean(lons), width=1.3*(max(lons) - min(lons))*1000*180, \
                         height=1.3*(max(lats) - min(lats))*1000*180, coords='geo',datetime=times[0])
    # (the 1000* factors are to replace 1000 in how usually these are written as "width=111e3*180")

    x,y = m(lons,lats,coords='geo')
    m.plot(x,y,'b',label="ePOP ground track")

    my_xticks = []
    num_ticks = 5
    length = times.__len__()
    tick_sep = length/(num_ticks - 1)
    alt_t = alts[0]
    lon_t = geog_longs[0]
    lat_t = geog_lats[0]
    dt_t  = times[0]
    my_xticks.append("Altitude:    "+str(alt_t)+"\nLatitude:    "+str(lat_t)+"\nLongitude:    "+str(lon_t)+"\nTime (UTC):    "+str(dt_t))
    for i in range(num_ticks-1):
        alt_t = alts[tick_sep*(i+1)]
        lon_t = geog_longs[tick_sep*(i+1)]
        lat_t = geog_lats[tick_sep*(i+1)]
        dt_t  = times[tick_sep*(i+1)]
        my_xticks.append(str(alt_t)+"\n"+str(lat_t)+"\n"+str(lon_t)+"\n"+str(dt_t))

    plt.xlabel('Geographic Longitude (degrees)')
    plt.ylabel('Geographic Latitude (degrees)')
    plt.title("EPOP Closest Approach vs. Ottawa radar for " + "2016-04-" + str(times[0].day))
    #plt.legend(loc='best')
    plt.show()
 

def plot_fov_sat(fovname, date, ephem_lons, ephem_lats, suppress_show=False):
    """
    This function uses matplotlib to conveniently plot a SuperDARN radar FOV 
    along with the EPOP's geographic coordinates, and saves the figure to an
    output directory (which is currently only checked for in script.py).

    ** ARGS **
        fovname (String): the name of the radar (as in DaVitPy, e.g. "Saskatoon")
        date (Datetime object): the desired date of the plot. Used for the title.
        ephem_lons (Numpy Array): Array of longitude points for a ground track.
        ephem_lats (Numpy Array): Array of latitude points for a ground track.
        [suppress_show] (Boolean): Tells function to just save the image as a file.
    
    ** RETURNS **
        - 
    """
    fov = get_fov_by_name(fovname)
    fovlons = ((fov.lonFull+360.)%360.).ravel()
    fovlats = (fov.latFull).ravel()
    fovcol = 3*np.ones(np.size(fovlons)) # 3 for blue FOVs

    # A workaround to list the blue FOV lines in the legend without listing 
    # them 17 separate times: give just one of them a label for the legend.
    lon = (fov.lonFull[0]+360.)%360.
    lat = (fov.latFull[0])
    plt.plot(lon,lat,'b',label="FOV")
    for i in range(np.shape(fov.lonFull)[0] - 1):
        lon = (fov.lonFull[i+1]+360.)%360.
        lat = (fov.latFull[i+1])
        f = plt.plot(lon,lat,'b')

    ephemcol = 5*np.ones(np.size(ephem_lats)) # 5 for a red track
    ephemlons = (ephem_lons+360.)%360.
    ephemlats = ephem_lats

    plt.title(fovname + " Radar FOV vs RRI Ephemeris for " + date.__str__())
    plt.xlabel('Geographic Longitude (degrees)')
    plt.ylabel('Geographic Latitude (degrees)')

    # So that the formats match, I will ensure months and days are padded for
    # the output figure's name as well.
    month = "0" + str(date.month) if str(date.month).__len__() == 1 else str(date.month)
    day = "0" + str(date.day) if str(date.day).__len__() == 1 else str(date.day)
    hr = "0" + str(date.hour) if str(date.hour).__len__() == 1 else str(date.hour)
    mn = "0" + str(date.minute) if str(date.minute).__len__() == 1 else str(date.minute)

    plt.plot(ephemlons,ephemlats,'r',label="RRI Ephemeris")
    plt.legend()
    plt.savefig("./data/output/"+str(date.year)+month+day+"_"+hr+"h"+mn+"_"+fovname) 
    if suppress_show==False:
        plt.show()
    
def get_fov_by_name(name):
    """
    This function shortcuts the FOV creation process by handling assumptions I
    typically make when creating FOV objects.

    ** ARGS **
        name (string): the name (as DaVitPy understands them) of a SuperDARN 
            radar whose FOV is to be returned.
    
    ** RETURNS **
        fov (pydarn.radar.radFov.fov): a (front) FOV object for the SuperDARN
            site with the given name.
    """
    from davitpy import pydarn
    nw = pydarn.radar.network()
    rad = nw.getRadarByName(name)
    assert (rad != False) #if name is illegal, getRadarByName returns False
    site = pydarn.radar.site(radId = rad.id)
    fov = pydarn.radar.radFov.fov(site=site,altitude=300.0,model='IS',coords='geo',ngates=75)
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

    
    plot_sat_ephemeris(date_string="20160418")
