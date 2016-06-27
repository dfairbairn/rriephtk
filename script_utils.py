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
from datetime import datetime

def ephems_to_datetime(ephem_times):
    """
    This function allows a whole array of ephemeris times to be conveniently
    converted in the same manner as with the ephem_to_datetime() function
    """
    assert isinstance(ephem_times, np.ndarray), "Not an array"
 
    # i) Check the seconds between May 24 1968 (ephem MET) and Jan 1 1970 (neg number)
    t_off = subprocess.check_output(["date", "--date=1968-05-24 0:00:00", "+%s"])
    t_off = float(t_off.split("\n",1)[0]) # extract the integer value

    # ii) Do the math.
    ephtimes = ephem_times + t_off
    times = []
    for i in range(np.size(ephtimes)):
        times.append(datetime.utcfromtimestamp(ephtimes[i]))
    return times

def ephem_to_datetime(ephem):
    """
    This function exists in order to conveniently convert the weird ephemeris
    time data format from the EPOP RRI instrument into a datetime object.

    *note*: makes use of the bash command 'date', so make sure you have it.

    """    
    # i) Check # seconds between May 24 1968 (ephem MET) and Jan 1 1970 (neg number)
    t_off = subprocess.check_output(["date", "--date=1968-05-24 0:00:00", "+%s"])
    t_off = float(t_off.split("\n",1)[0]) # extract the integer value
  
    # ii) Do the math. 
    return datetime.utcfromtimestamp(t_off + float(ephem))    


#TODO: plot the back lobe?
def plot_fov_sat(fovname, date, ephem_lons, ephem_lats,suppress_show=False):
    """
    This function uses matplotlib to conveniently plot a SuperDARN radar FOV 
    along with the EPOP's geographic coordinates, and saves the figure to CWD.

    **PARAMS**
    fovname: (String) the name of the radar (as in DaVitPy, e.g. "Saskatoon")
    date: (Datetime object) the desired date of the plot. Used for the title.
    ephem_lons: (Numpy Array) Array of longitude points for a ground track.
    ephem_lats: (Numpy Array) Array of latitude points for a ground track.
    [suppress_show]: (Boolean) Tells function to just save the image as a file.

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

    plt.plot(ephemlons,ephemlats,'r',label="RRI Ephemeris")
    plt.legend()
    plt.savefig("./data/output/"+str(date.year)+"_"+month+day+"_"+fovname) 
    if suppress_show==False:
        plt.show()
    


def plot_fovs_sat(fovnames, date, ephem_longs, ephem_lats):
    """
    This function exists in order to facilitate plotting multiple FOVs together
    with the ephemeris data.
    #TODO: update like with the above function
    """

    if (np.shape(fovnamess)[0] > 10):
        print "Can't do more than 10 FOVs"
        return


    lons = []
    lats = []
    colors = []
    color_offset = 0
    for fovname in fovnames:
        fov = get_fov_by_name(fovname)
        fovlons = ((fov.lonFull+360.)%360.).ravel()
        fovlats = (fov.latFull).ravel()
        fovcol = (1+color_offset)*np.ones(np.size(fovlons))
        color_offset += 1
        
        lons = np.concatenate((lons,fovlons))
        lats = np.concatenate((lats,fovlats))
        colors = np.concatenate((colors,fovcol))
    
    ephemlons = (ephem_longs+360.)%360.
    ephemlats = ephem_lats 
    ephemcol = (color_offset+1)*np.ones(np.size(ephemlons))
    lons = np.concatenate((lons,ephemlons)) 
    lats = np.concatenate((lats,ephemlats))
    colors = np.concatenate((colors,ephemcol))
 
    plt.title("Radar FOVs vs RRI Ephemeris for " + date.__str__())
    plt.xlabel('Geographic Longitude (degrees)')
    plt.ylabel('Geographic Latitude (degrees)')
   
    plt.plot(lons,lats,c=colors)#,edgecolors='face')
    plt.show()


def get_fov_by_name(name):
    """
    This function shortcuts the FOV creation process by handling assumptions I
    typically make when creating FOV objects.
    """
    from davitpy import pydarn
    nw = pydarn.radar.network()
    rad_id = (nw.getRadarByName(name)).id
    site = pydarn.radar.site(radId = rad_id)
    fov = pydarn.radar.radFov.fov(site=site,altitude=300.0,model='IS',coords='geo',ngates=75)
    return fov

def two_pad(in_time):
    """ 
    Takes in a number of 1 or 2 digits, returns a string of two digits. 
    """
    return "0" + str(in_time) if str(in_time).__len__() == 1 else str(in_time)

