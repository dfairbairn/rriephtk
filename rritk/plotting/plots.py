"""
file: plots.py
description:
    A file which contains a variety of plotting procedures specific to
    analysis of RRI experiments (e.g. CASSIOPE ground tracks, aspect
    angle from transmitter to receiver, etc.).
author: David Fairbairn
date: July 2017

"""

import subprocess 
import datetime as dt #A wrapper class for 'file' which tracks line numbers. 

import numpy as np 
import matplotlib.pyplot as plt

from davitpy.utils import plotUtils

import __init__
import rritk.utils.data_utils as datautils

OTTAWA_TX_LON = -75.552
OTTAWA_TX_LAT = 45.403
OTTAWA_TX_ELEV = 0.070 # 70m elevation in Ottawa - small compared to 300km satellite
MILLSTONE_TX_LON = -71.491 # acquired from google maps satellite imagery
MILLSTONE_TX_LAT = 42.619
data_path = '../rri-conjunction-script/data'

def plot_sat_ephemeris(date_string=None,lons=None,lats=None,alts=None,ephtimes=None):
    """


    """
    
    if lons==None or lats==None or alts==None or ephtimes==None:
        if date_string==None or not isinstance(date_string, type("e.g.")):
            print "Need to provide either lons/lats/alts/ephtimes or date_string!"
            return -1
        # e.g. like paired with the MGF data access
        rri_fname,indx_rev = datautils.get_ottawa_data(date_string) 
        lons,lats,alts,ephtimes = datautils.get_rri_ephemeris(rri_fname)
    times=  datautils.ephems_to_datetime(ephtimes)
    
    # A different font for the legend etc. might be nice
    #fig = plt.figure()
    font = {'fontname':'Computer Modern'}
    m = plotUtils.mapObj(lat_0=np.mean(lats), lon_0=np.mean(lons), width=4.0*(max(lons) - min(lons))*1000*180, \
                         height=1.3*(max(lats) - min(lats))*1000*180, coords='geo',resolution='i',datetime=times[0])
    # (the 1000* factors are to replace 1000 in how usually these are written as "width=111e3*180")

    x,y = m(lons,lats,coords='geo')
    m.plot(x,y,'b',label="ePOP ground track")

    my_xticks = []
    num_ticks = 5
    length = times.__len__()
    tick_sep = length/(num_ticks - 1)
    alt_t = alts[0]
    lon_t = lons[0]
    lat_t = lats[0]
    dt_t  = times[0]
    my_xticks.append("Altitude:    "+str(alt_t)+"\nLatitude:    "+str(lat_t)+"\nLongitude:    "+str(lon_t)+"\nTime (UTC):    "+str(dt_t))
    for i in range(num_ticks-1):
        alt_t = alts[tick_sep*(i+1)]
        lon_t = lons[tick_sep*(i+1)]
        lat_t = lats[tick_sep*(i+1)]
        dt_t  = times[tick_sep*(i+1)]
        my_xticks.append(str(alt_t)+"\n"+str(lat_t)+"\n"+str(lon_t)+"\n"+str(dt_t))

    plt.xlabel('Geographic Longitude (degrees)')
    plt.ylabel('Geographic Latitude (degrees)')
    plt.title("EPOP Closest Approach vs. Ottawa radar for " + "2016-04-" + str(times[0].day))
    #plt.legend(loc='best')
    plt.show()
 

def plot_fov_sat(fovname, ephem_lons, ephem_lats, date=None, frontback='front',
                 suppress_show=False, outfile_name=None):
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
    fov = datautils.get_fov_by_name(fovname, frontback=frontback)
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

    t = fovname + " Radar FOV vs RRI Ephemeris" 
    if date is not None:
        t = t + " for " + str(date)
        # So that the formats match, I will ensure months and days are padded for
        # the output figure's name as well.
        month = "0" + str(date.month) if str(date.month).__len__() == 1 else str(date.month)
        day = "0" + str(date.day) if str(date.day).__len__() == 1 else str(date.day)
        hr = "0" + str(date.hour) if str(date.hour).__len__() == 1 else str(date.hour)
        mn = "0" + str(date.minute) if str(date.minute).__len__() == 1 else str(date.minute)


    plt.title(t)
    plt.xlabel('Geographic Longitude (degrees)')
    plt.ylabel('Geographic Latitude (degrees)')

    plt.plot(ephemlons,ephemlats,'r',label="RRI Ephemeris")
    plt.legend()
    if outfile_name is not None:
        plt.savefig(outfile_name)
    elif date is not None:
        plt.savefig("./data/output/"+str(date.year)+month+day+"_"+hr+"h"+mn+"_"+fovname+".png") 
    else:
        plt.savefig("someplot.png")
    if suppress_show==False:
        plt.show()

def ephem_ticks(lons,lats,alts,ephtimes,mlons,mlats,mlts):
    """
    When you're in the middle of plotting some RRI data for which the data is once 
    a second like the basic ephemeris data, this prepares the nice x-axis with all ephem
    info (so you just continue saying 'plt.plot(x,y), plt.show()' after calling this)
    """
    times = ephems_to_datetime(ephtimes)
    my_xticks = []
    num_ticks = 5
    length = times.__len__()
    tick_sep = length/(num_ticks - 1)
    dt_t  = times[0]
    alt_t = alts[0]
    lon_t = lons[0]
    lat_t = lats[0]
    mlon_t = mlons[0]
    mlat_t = mlats[0]
    mlt_t = mlts[0]
    my_xticks.append("Time (UTC):    "+str(dt_t.time())+"\nLatitude:    "+str(lat_t)+\
        "\nLongitude:    "+str(lon_t)+"\nAltitude:    "+str(alt_t)+\
    "\nMagnetic Local Time:    "+str(mlt_t)+"\nMagnetic Latitude:    "+str(mlat_t)+\
    "\nMagnetic Longitude:    "+str(mlon_t))
    for i in range(num_ticks-1):
        alt_t = alts[tick_sep*(i+1)]
        lon_t = lons[tick_sep*(i+1)]
        lat_t = lats[tick_sep*(i+1)]
        mlon_t = mlons[tick_sep*(i+1)]
        mlat_t = mlats[tick_sep*(i+1)]
        dt_t  = times[tick_sep*(i+1)]
        mlt_t = mlts[tick_sep*(i+1)]
        my_xticks.append(str(dt_t.time())+"\n"+str(lat_t)+"\n"+str(lon_t)+"\n"+str(alt_t)+\
                "\n"+str(mlt_t)+"\n"+str(mlat_t)+"\n"+str(mlon_t))
    
    indices = range(times.__len__())
    tick_indices = [i*tick_sep for i in range(num_ticks)]
    plt.xticks(tick_indices, my_xticks)
   
if __name__=="__main__":
    path, fname = datautils.initialize_data()
    lons, lats, alts, ephtimes = datautils.get_rri_ephemeris(fname)
    plot_fov_sat("Saskatoon", lons, lats, date=datautils.ephem_to_datetime(ephtimes[0])) 
