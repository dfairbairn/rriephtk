"""
file: 'ottawa_plots.py'
author: David Fairbairn
date: July 2016
description: This file contains code used to produce fly-by plots of CASSIOPE
    vs the Ottawa-based radar that was running an RRI experiment on the dates
    April 18-22nd.

    On each day, faraday rotation behaviour was seen to be occurring, reaching
    a point where it would stop rotating in one direction and start rotating
    the opposite way. This turned out to _NOT_ be a simple case of the distance
    from CASSIOPE to Ottawa reaching a minimum at this same point, so these 
    plots serve to contrast the location of the rotation inversion, and the 
    point of closest approach.

    In addition, magnetic coordinates and IGRF can be included.

"""
from script_utils import *
from data_utils import *

from davitpy.utils import plotUtils
from davitpy.models import aacgm

from datetime import datetime as dt
from davitpy.models import igrf

import math
#import matplotlib.pyplot as plt
#import numpy as np
import sys

def get_ottawa_kvec(glon, glat, altitude, time):
    """
    This function takes a satellite ephemeris point as input, calculating
    the straight-line k-vector from the Ottawa-based transmitter.
    """
    # The specific coordinates of the NRCAN geomagnetic    
    ottawa_lon = -75.552 
    ottawa_lat = 45.403

    kx,ky,kz = (0.,0.,0.)
    return (kx,ky,kz)

def get_bvec(glon, glat, altitude, time):
    """
    This function uses the IGRF model in DavitPy to calculate the magnetic
    field vector at a given near-Earth location at a particular time.

    Currently, the function just takes geographic longitude and latitude. It 
    shouldn't be difficult to extend its functionality to accept magnetic
    coordinates as well thoug.

    **PARAMS**
    glon (Float): Geographic longitude of point
    glat (Float): Geographic latitude of point
    alt (Float): Altitude from Earth's surface in km
    time (Datetime): Time of interest (magnetic North's location varies year-to-year)
    """
    from davitpy.models import igrf
    from davitpy import utils
    itype = 1 #Geodetic coordinates
    date = utils.dateToDecYear(time) # decimal year
    alt = altitude
    stp = 1. #
    # The IGRF function takes a grid of latitude and longitude points for which
    # to calculate the field. We just want one point.
    xlti, xltf, xltd = glat, glat, stp
    xlni, xlnf, xlnd = glon, glon, stp 
    ifl = 0 # Main field
    # Call fortran subroutine
    lat,lon,d,s,h,bx,by,bz,f = igrf.igrf11(itype,date,alt,ifl,xlti,xltf,xltd,xlni,xlnf,xlnd)
    return (bx,by,bz)

def plot_ottawa_ephem(date_string):
    """
    Put the plotting procedure for looking at satellite ephemeris vs. Ottawa
    transmitter into a function.

    Note that the transmitter is a Barker & Williamson Model 110.

    **PARAMS**
    date_string (String): String in the format of "20160418"

    """
    if isinstance(date_string, NoneType): date_string="20160418"
       
    # **** CHOOSE ONE OF THESE RRI FILES THEN RUN THE SCRIPT ****
    if "20160418"==date_string:
        geog_longs,geog_lats,alts,ephemtimes = get_rri_ephemeris("./data/RRI_20160418_222759_223156_lv1_v2.h5") #18th
        index_inversion = 167 #for 18th
    elif "20160419"==date_string:
        index_inversion = 178 #for 19th
        geog_longs,geog_lats,alts,ephemtimes = get_rri_ephemeris("./data/RRI_20160419_220939_221336_lv1_v2.h5") #19th
    elif "20160420"==date_string:
        index_inversion = 213 #for 20th
        geog_longs,geog_lats,alts,ephemtimes = get_rri_ephemeris("./data/RRI_20160420_215117_215514_lv1_v2.h5") #20th
    elif "20160421"==date_string:
        index_inversion = 205 #?? for 21st?
        geog_longs,geog_lats,alts,ephemtimes = get_rri_ephemeris("./data/RRI_20160421_213255_213652_lv1_v2.h5") #21st
    elif "20160422"==date_string:
        index_inversion = 222 #for 22nd
        geog_longs,geog_lats,alts,ephemtimes = get_rri_ephemeris("./data/RRI_20160422_211435_211832_lv1_v2.h5") #22nd
    else:
        print "Invalid input date."
        return None    

    # Location of Ottawa: I looked it up and am hard-coding it here.
    ottawa_long = -75.552
    ottawa_lat = 45.403 
    times = ephems_to_datetime(ephemtimes)

    # Make all longitudes positive?
    #ottawa_long = (ottawa_long+360.)%360.0
    #geog_longs = (geog_longs+360.)%360.0

    # *** FINDING THE CLOSEST APPROACH ***
    # Using the Haversine formula (in a function in script_utils.py), the closest
    # approach is determined by brute force.
    dists = []
    longdists = []
    latdists = []
    shortest_dist = sys.maxint #Initially set this very high
    for i in range(np.size(geog_longs)):
        dist = haversine(geog_longs[i], geog_lats[i], ottawa_long, ottawa_lat)
        # Initially I took a quick and dirty approach to find the point of smallest
        # Euclidean distance in terms of latitudes and longitudes.
        #longdist = abs(geog_longs[i] - ottawa_long) # difference of longitudes
        #latdist = abs(geog_lats[i] - ottawa_lat) # difference of latitudes
        #dist = np.sqrt(longdist*longdist + latdist*latdist)  
        if dist < shortest_dist:
            shortest_dist = dist
            shortest = i

    appr_time = times[shortest]
    
    # The numeric data type that I was retrieving from geog_longs, when _NOT_ stored
    # in an array, was being rejected by the mapObj() function below. So I convert 
    # these numbers to floats explicitly here.
    shortest_ephem_long = float(geog_longs[shortest])
    shortest_ephem_lat = float(geog_lats[shortest])
    inversion_ephem_long = float(geog_longs[index_inversion])
    inversion_ephem_lat = float(geog_lats[index_inversion])
    
    # A different font for the legend etc. might be nice
    #fig = plt.figure()
    font = {'fontname':'Computer Modern'}
    m = plotUtils.mapObj(lat_0=45.0, lon_0=-75.0, width=111e3*180, height=111e3*90, coords='geo',datetime=times[0])
    
    # FIRST: Plot the location of Ottawa
    x,y = m(ottawa_long,ottawa_lat,coords='geo')
    m.plot(x,y,'ro',label="Ottawa")
    
    # SECOND: Plot the satellite ground-track.
    x,y = m(geog_longs, geog_lats, coords='geo')
    m.plot(x,y,'b',label="EPOP ground track")
    
    # THIRD: Plot a circle emphasizing the point of closest approach
    x,y = m(shortest_ephem_long,shortest_ephem_lat, coords='geo')
    m.plot(x,y,'bo',label=("Shortest Approach at " + str(appr_time)))
    
    # FOURTH: Plot the line from Ottawa to the nearest approach of the satellite.
    x,y = m([shortest_ephem_long, ottawa_long], [shortest_ephem_lat, ottawa_lat], coords='geo')
    m.plot(x,y,'g')
    
    # FIFTH: Plot the piont I've determined is the point of the Faraday Rotation inversion.
    x,y = m(inversion_ephem_long,inversion_ephem_lat,coords='geo')
    m.plot(x,y,'yo',label=("Inversion of Faraday Rotation"))
    
    # SIXTH: a few lines of magnetic longitude and latitude will be plotted as well.
    """
    merid1_mlat = merid2_mlat = merid3_mlat = np.arange(181) - 90.
    merid1_mlon = merid1_mlat*0. - 10.
    merid2_mlon = merid2_mlat*0.
    merid3_mlon = merid3_mlat*0. + 10.
    paral1_mlon = paral2_mlon = paral3_mlon = np.arange(361) - 45.
    paral1_mlat = paral1_mlon*0. + 55.
    paral2_mlat = paral2_mlon*0. + 60.
    paral3_mlat = paral3_mlon*0. + 65.
    
    zero_alts = merid2_mlon
    
    merid1_glat,merid1_glon,r = aacgm.aacgmConvArr(merid1_mlat.tolist(),merid1_mlon.tolist(),zero_alts.tolist(),2016,1)
    merid2_glat,merid2_glon,r = aacgm.aacgmConvArr(merid2_mlat.tolist(),merid2_mlon.tolist(),zero_alts.tolist(),2016,1)
    merid3_glat,merid3_glon,r = aacgm.aacgmConvArr(merid3_mlat.tolist(),merid3_mlon.tolist(),zero_alts.tolist(),2016,1)
    paral1_glat,paral1_glon,r = aacgm.aacgmConvArr(paral1_mlat.tolist(),paral1_mlon.tolist(),zero_alts.tolist(),2016,1)
    paral2_glat,paral2_glon,r = aacgm.aacgmConvArr(paral2_mlat.tolist(),paral2_mlon.tolist(),zero_alts.tolist(),2016,1)
    paral3_glat,paral3_glon,r = aacgm.aacgmConvArr(paral3_mlat.tolist(),paral3_mlon.tolist(),zero_alts.tolist(),2016,1)
    
    x,y = m(merid1_glon, merid1_glat, coords='mag')#coords='geo')
    m.plot(x,y,'k')#,label="Line of Magnetic Longitude of -20 Degrees")
    
    x,y = m(merid2_glon, merid2_glat, coords='mag')#coords='geo')
    m.plot(x,y,'k')#,label="Line of Magnetic Longitude of 0 Degrees")
    
    x,y = m(merid3_glon, merid3_glat, coords='mag')#coords='geo')
    m.plot(x,y,'k')#,label="Line of Magnetic Longitude of +20 Degrees")
    
    x,y = m(paral1_glon, paral1_glat, coords='mag')#coords='geo')
    m.plot(x,y,'k')#,label="Line of Magnetic Latitude of +35 Degrees")
    
    x,y = m(paral2_glon, paral2_glat, coords='mag')#coords='geo')
    m.plot(x,y,'k')#,label="Line of Magnetic Latitude of +55 Degrees")
    
    x,y = m(paral3_glon, paral3_glat, coords='mag')# coords='geo')
    m.plot(x,y,'k')#,label="Line of Magnetic Latitude of +75 Degrees")
    """
    
    N = 10
    latdivs = 90./N
    londivs = 360./N
    for n in range(N):
        merid_mlat = np.arange(181) - 90.
        merid_mlon = merid_mlat*0 - n*londivs
        paral_mlon = np.arange(361) - 180.
        paral_mlat = paral_mlon*0. + n*latdivs
        x,y = m(merid_mlon, merid_mlat, coords='mag')
        m.plot(x,y,'k')#,label="Line of Magnetic Longitude of -20 Degrees")
        x,y = m(paral_mlon, paral_mlat, coords='mag')
        m.plot(x,y,'k')#,label="Line of Magnetic Latitude of +75 Degrees")
    
    # SEVENTH: GET IGRF DATA FOR EACH EPHEMERIS POINT
    """
    itype = 1 #Geodetic coordinates
    pyDate = times[0] # The first time we pull from the RRI file
    date = utils.dateToDecYear(pyDate) # decimal year
    alt = 300. # altitude #TODO: grab altitudes for series of satellite positions we care about.
    stp = 1. #
    xlti, xltf, xltd = ottawa_lat, ottawa_lat,stp # latitude start, stop, step
    xlni, xlnf, xlnd = ottawa_long, ottawa_long,stp # longitude start, stop, step
    ifl = 0 # Main field
    # Call fortran subroutine
    lat,lon,d,s,h,x,y,z,f = igrf.igrf11(itype,date,alt,ifl,xlti,xltf,xltd,xlni,xlnf,xlnd)
    """

    plt.xlabel('Geographic Longitude (degrees)')
    plt.ylabel('Geographic Latitude (degrees)')
    plt.title("EPOP Closest Approach vs. Ottawa radar for " + "2016-04-" + str(times[0].day))
    #plt.legend(loc='best')
    plt.show()

# -----------------------------------------------------------------------------
ottawa_long = -75.552
ottawa_lat = 45.403
lon = ottawa_long
lat = ottawa_lat
alt = 0.
(bx,by,bz) = get_bvec(lon,lat,alt,dt.now())

