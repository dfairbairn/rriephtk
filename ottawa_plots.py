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
import matplotlib.pyplot as plt
import numpy as np
import sys

OTTAWA_TX_LON = -75.552
OTTAWA_TX_LAT = 45.403
OTTAWA_TX_ELEV = 0.070 # 70m elevation in Ottawa - small compared to 300km satellite
def get_ottawa_kvec(glon, glat, altitude, time):
    """
    This function takes a satellite ephemeris point as input, calculating
    the straight-line k-vector from the Ottawa-based transmitter in terms of 
    North (x), East (y), Down (z) components.
    """
    #TODO: Parallelize/accept vector inputs

    # The specific coordinates of the NRCAN geomagnetic    
    OTTAWA_TX_LON = -75.552 
    OTTAWA_TX_LAT = 45.403

    # In spherical coordinates, subtracting vectors doesn't get us the path
    # from one point to another along the surface. To accomplish that, we use
    # bearings and elevation angles.
    
    init_bearing,final_bearing = get_bearing(OTTAWA_TX_LON, OTTAWA_TX_LAT, glon, glat)
    elev_angle = get_elevation_angle(OTTAWA_TX_LON, OTTAWA_TX_LAT, OTTAWA_TX_ELEV, glon, glat, altitude)

    # for calculation of vector components, we need theta (the bearing in the
    # N-E plane) at CASSIOPE, and phi, the angle-from-down (toward the N-E plane)
    theta = final_bearing
    phi = 90. + elev_angle
    kx = np.sin(np.deg2rad(phi))*np.cos(np.deg2rad(theta))
    ky = np.sin(np.deg2rad(phi))*np.sin(np.deg2rad(theta))
    kz = np.cos(np.deg2rad(phi)) 
    return np.array((kx,ky,kz))

def get_bearing(lon1, lat1, lon2, lat2):
    """
    Gives the bearing clockwise from north in degrees from point 1 to point2,
    at point 1 (init_bearing) and at point 2 (final_bearing)
    """
    #TODO: Parallelize/accept vector inputs
    atan2_argx = np.sin(np.deg2rad(lon2 - lon1))*np.cos(np.deg2rad(lat2))
    atan2_argy1 = np.cos(np.deg2rad(lat1))*np.sin(np.deg2rad(lat2))
    atan2_argy2 = np.sin(np.deg2rad(lat1))*np.cos(np.deg2rad(lat2))*np.cos(np.deg2rad(lon2 - lon1))
    init_bearing = np.rad2deg(np.arctan2(atan2_argx, atan2_argy1 - atan2_argy2))

    # The FINAL bearing is just the bearing from final point to the initial 
    # point, reversed by adding 180 modulo 360
    atan2_argx = np.sin(np.deg2rad(lon1 - lon2))*np.cos(np.deg2rad(lat1))
    atan2_argy1 = np.cos(np.deg2rad(lat2))*np.sin(np.deg2rad(lat1))
    atan2_argy2 = np.sin(np.deg2rad(lat2))*np.cos(np.deg2rad(lat1))*np.cos(np.deg2rad(lon1 - lon2))
    rev_bearing = np.rad2deg(np.arctan2(atan2_argx, atan2_argy1 - atan2_argy2))
    final_bearing = (rev_bearing + 180.) % 360.
    return init_bearing,final_bearing

def get_elevation_angle(lon1, lat1, alt1, lon2, lat2, alt2):
    """
    
    """
    arcdist = haversine(lon1, lat1, lon2, lat2)
    delta_alt = alt2 - alt1
    elev_angle = np.rad2deg(np.arctan2(delta_alt,arcdist))
    return elev_angle

def get_bvec(glon, glat, altitude, time):
    """
    This function uses the IGRF model in DavitPy to calculate the magnetic
    field vector at a given near-Earth location at a particular time.

    Currently, the function just takes geographic longitude and latitude. It 
    shouldn't be difficult to extend its functionality to accept magnetic
    coordinates as well though.

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
    return np.array((bx[0],by[0],bz[0]))

def get_kb_ottawa_angle(ephem_glongs, ephem_glats, ephem_alts, ephem_times):
    """
    This function takes the ephemeris data as arguments and using the Ottawa
    NRCAN transmitter, compares the k (line of sight) vector and the IGRF B 
    field vector to determine the relative angle between the propagating radio
    wave and the B-field in the ionospheric plasma.

    The purpose of this is to get an idea of which mode the radio wave is 
    propagating under (and if it's undergoing Faraday rotation, etc.)

    *** PARAMS ***
    ephem_glongs (np.array of floats): ground-track geographic longitude
    ephem_glats (np.array of floats): ground-track geographic latitude
    ephem_alts (np.array of floats): altitude of CASSIOPE 
    ephem_times (np.array of floats): in truncated JD (MET), seconds since era

    *** RETURNS ***
    bvecs (np.array of float triplets): the IGRF B field at each ephemeris point
    kvecs (np.array of float triplets): the tx-to-rx(sat) vector at each ephemeris point
    angles (np.array of floats): the aspect angle (angle between respective bvecs and kvecs)
    """
    times = ephems_to_datetime(ephem_times)
    angles = []
    bvecs = []
    kvecs = []
    for i in range(ephem_glongs.__len__()):
        lon = ephem_glongs[i]
        lat = ephem_glats[i]
        alt = ephem_alts[i]
        time = times[i]
        bvec = get_bvec(lon,lat,alt,time)
        kvec = get_ottawa_kvec(lon,lat,alt,time)
        #print "B vector: " + str(bvec)
        #print "K vector: " + str(kvec)
        # Take the dot product, divide out the magnitude of the vectors
        prod = np.dot(kvec,bvec)/(np.linalg.norm(bvec)*np.linalg.norm(kvec))
        angle = np.rad2deg(np.arccos(prod)) 
        angles.append(angle)
        bvecs.append(bvec)
        kvecs.append(kvec) 
    return bvecs,kvecs,angles

def get_ramdirs(glon, glat, altitude, time):
    """
    Computes the velocity components of the satellite based on its ephemeris.
    
    *** PARAMS ***
    glon (np.array of floats): the longitudes of the satellite (deg) 
    glat (np.array of floats): the latitudes of the satellite (deg)
    altitude (np.array of floats): the altitudes of the satellite (deg)
    time ( doesnt even matter right now ) : The times (currently unused)

    *** RETURNS ***
    v (list of floats): velocity vector components (in km/s)
    dists (list of floats): full distances traveled in 1 sec in km

    *** TESTING? ***
    calculated velocity values agreed with velocity values from other sources


    """
    v = []
    dists = [] 
    for i in range(glon.__len__()-1):
        lon1 = glon[i]
        lon2 = glon[i+1]
        lat1 = glat[i]
        lat2 = glat[i+1]
        alt1 = altitude[i]
        alt2 = altitude[i+1]
        latdist = haversine(lon1,lat1,lon1,lat2) # if lon1,lat1 -- lon2,lat2 is the hypotenuse, this is the vertical
        londist = haversine(lon1,lat1,lon2,lat1) # and this is the horizontal
        altdist = alt2 - alt1
        dist = np.linalg.norm((latdist,londist,altdist))
        v.append((latdist,londist,altdist))
        dists.append(dist)
    return v,dists

def get_closest_ottawa_approach(glons, glats, alts):
    """
    Uses a haversine formula-based distance calculation to find the closest 
    approach of a set of points describing a ground-track.

    *** PARAMS *** 
    glons (np.array of floats): longitudes (deg) of ground-track points
    glats (np.array of floats): latitudes (deg) of ground-track points
    alts (np.array of floats): altitudes (in km)

    ** RETURNS **
    indx_shortest (integer): index of the lat/lon pair where the ottawa transmitter is closest
    dist_shortest (float): the actual closest distance
    """

    # *** FINDING THE CLOSEST APPROACH ***
    # Using the Haversine formula (in a function in script_utils.py), the closest
    # approach is determined by brute force.
    dists = []
    longdists = []
    latdists = []
    dist_shortest = sys.maxint #Initially set this very high
    for i in range(np.size(glons)):
        dist_before_alt = haversine(glons[i], glats[i], OTTAWA_TX_LON, OTTAWA_TX_LAT)
        delt_alt = alts[i] - OTTAWA_TX_ELEV
        dist = np.linalg.norm((dist_before_alt,delt_alt))
        # Initially I took a quick and dirty approach to find the point of smallest
        # Euclidean distance in terms of latitudes and longitudes.
        #longdist = abs(geog_longs[i] - OTTAWA_TX_LON) # difference of longitudes
        #latdist = abs(geog_lats[i] - OTTAWA_TX_LAT) # difference of latitudes
        #dist = np.sqrt(longdist*longdist + latdist*latdist)  
        dists.append(dist)
        if dist < dist_shortest:
            dist_shortest = dist
            indx_shortest = i
    return indx_shortest, dists 

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
        fname = "./data/RRI_20160418_222759_223156_lv1_v2.h5" #18th
        index_reversal = 167 #for 18th # Stay
    elif "20160419"==date_string:
        #index_reversal = 178 #for 19th # Could go -1
        index_reversal = 177
        fname = "./data/RRI_20160419_220939_221336_lv1_v2.h5" #19th
    elif "20160420"==date_string:
        index_reversal = 213 #for 20th # Stay
        fname = "./data/RRI_20160420_215117_215514_lv1_v2.h5" #20th
    elif "20160421"==date_string:
        #index_reversal = 205 #?? for 21st? # Could go up +5
        index_reversal = 210
        fname = "./data/RRI_20160421_213255_213652_lv1_v2.h5" #21st
    elif "20160422"==date_string:
        #index_reversal = 222 #for 22nd # Could go down -1 
        index_reversal = 221
        fname = "./data/RRI_20160422_211435_211832_lv1_v2.h5" #22nd
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
        fname = "./data/RRI_20160418_222759_223156_lv1_v2.h5" #18th
        orientation_rev = 168 #for 18th # Stay
        ellipt_rev = 167
    elif "20160419"==date_string:
        orientation_rev = 177
        ellipt_rev = 158
        fname = "./data/RRI_20160419_220939_221336_lv1_v2.h5" #19th
    elif "20160420"==date_string:
        orientation_rev  = 213 #for 20th # Stay
        ellipt_rev = 133 # ??? it also does something near 213 though
        fname = "./data/RRI_20160420_215117_215514_lv1_v2.h5" #20th
    elif "20160421"==date_string:
        orientation_rev  = 210
        ellipt_rev = 125
        fname = "./data/RRI_20160421_213255_213652_lv1_v2.h5" #21st
    elif "20160422"==date_string:
        orientation_rev  = 221
        ellipt_rev = 110 
        fname = "./data/RRI_20160422_211435_211832_lv1_v2.h5" #22nd
    else:
        print("Invalid input date.")
        return None    
    return fname,orientation_rev,ellipt_rev

   
 
def plot_ottawa_ephem(date_string):
    """
    Put the plotting procedure for looking at satellite ephemeris vs. Ottawa
    transmitter into a function.

    Note that the transmitter is a Barker & Williamson Model 110.

    **PARAMS**
    geog_longs
    geog_lats
    alts
    ephemtimes
    date_string (String): String in the format of "20160418"

    *** RETURNS ***
    - (Just plots)

    """
    # TODO: fixup this documentation
    if isinstance(date_string, type(None)): date_string="20160418"
       
    fname,index_reversal = get_ottawa_data(date_string)
    geog_longs,geog_lats,alts,ephemtimes = get_rri_ephemeris(fname)

    # Location of Ottawa: I looked it up and hard-coded it at the top
    times = ephems_to_datetime(ephemtimes)

    indx_shortest, dists = get_closest_ottawa_approach(geog_longs, geog_lats, alts)
    appr_time = times[indx_shortest]
    
    # The numeric data type that I was retrieving from geog_longs, when _NOT_ stored
    # in an array, was being rejected by the mapObj() function below. So I convert 
    # these numbers to floats explicitly here.
    shortest_ephem_long = float(geog_longs[indx_shortest])
    shortest_ephem_lat = float(geog_lats[indx_shortest])
    inversion_ephem_long = float(geog_longs[index_reversal])
    inversion_ephem_lat = float(geog_lats[index_reversal])
    
    # A different font for the legend etc. might be nice
    #fig = plt.figure()
    font = {'fontname':'Computer Modern'}
    m = plotUtils.mapObj(lat_0=45.0, lon_0=-75.0, width=111e3*180, height=111e3*90, coords='geo',datetime=times[0])
    
    # FIRST: Plot the location of Ottawa
    x,y = m(OTTAWA_TX_LON,OTTAWA_TX_LAT,coords='geo')
    m.plot(x,y,'r-o',markersize=9,label="Ottawa")
    
    # SECOND: Plot the satellite ground-track.
    x,y = m(geog_longs, geog_lats, coords='geo')
    m.plot(x,y,'b',label="EPOP ground track")
    
    # THIRD: Plot a circle emphasizing the point of closest approach
    x,y = m(shortest_ephem_long,shortest_ephem_lat, coords='geo')
    m.plot(x,y,'bo',label=("Closest Approach at " + str(appr_time)))
    
    # FOURTH: Plot the line from Ottawa to the nearest approach of the satellite.
    x,y = m([shortest_ephem_long, OTTAWA_TX_LON], [shortest_ephem_lat, OTTAWA_TX_LAT], coords='geo')
    m.plot(x,y,'g')
    
    # FIFTH: Plot the piont I've determined is the point of the Faraday Rotation inversion.
    x,y = m(inversion_ephem_long,inversion_ephem_lat,coords='geo')
    m.plot(x,y,'yo',label=("Inflection of Faraday Rotation"))
    
    # SIXTH: a few lines of magnetic longitude and latitude will be plotted as well.
   
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
    xlti, xltf, xltd = OTTAWA_TX_LAT, OTTAWA_TX_LAT,stp # latitude start, stop, step
    xlni, xlnf, xlnd = OTTAWA_TX_LON, OTTAWA_TX_LON,stp # longitude start, stop, step
    ifl = 0 # Main field
    # Call fortran subroutine
    lat,lon,d,s,h,x,y,z,f = igrf.igrf11(itype,date,alt,ifl,xlti,xltf,xltd,xlni,xlnf,xlnd)
    """

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

def plot_kb_angle(date_string):
    """
    Take date_string as input, produce kb-angle as output 

    *** PARAMS ***
    date_string (string): string in form "20160418" to denote date for which to plot KB angle

    *** RETURNS ***
    - (just plots)

    """
    fname,index_reversal = get_ottawa_data(date_string)
    lons,lats,alts,ephtimes,mlons,mlats,mlts,pitch,yaw,roll = get_rri_ephemeris_full(fname)
    bvecs,kvecs,angles = get_kb_ottawa_angle(lons,lats,alts,ephtimes)
    indx_closest, dists = get_closest_ottawa_approach(lons,lats,alts)
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
    dist_t = dists[0]
    my_xticks.append("Time (UTC):    "+str(dt_t.time())+"\nLatitude:    "+str(lat_t)+\
        "\nLongitude:    "+str(lon_t)+"\nAltitude:    "+str(alt_t)+\
    "\nMagnetic Local Time:    "+str(mlt_t)+"\nMagnetic Latitude:    "+str(mlat_t)+\
    "\nMagnetic Longitude:    "+str(mlon_t)+"\nDistance (km):    "+str(dist_t))
    for i in range(num_ticks-1):
        alt_t = alts[tick_sep*(i+1)]
        lon_t = lons[tick_sep*(i+1)]
        lat_t = lats[tick_sep*(i+1)]
        mlon_t = mlons[tick_sep*(i+1)]
        mlat_t = mlats[tick_sep*(i+1)]
        dt_t  = times[tick_sep*(i+1)]
        mlt_t = mlts[tick_sep*(i+1)]
        dist_t = dists[tick_sep*(i+1)]
        my_xticks.append(str(dt_t.time())+"\n"+str(lat_t)+"\n"+str(lon_t)+"\n"+str(alt_t)+\
                "\n"+str(mlt_t)+"\n"+str(mlat_t)+"\n"+str(mlon_t)+"\n"+str(dist_t))
    
    indices = range(angles.__len__())
    tick_indices = [i*tick_sep for i in range(num_ticks)]
    
    # Formatting adjustment so as to make room for the lengthy Ephemeris info
    fig = plt.figure()
    ax = plt.subplot(111,aspect = 'equal')
    plt.subplots_adjust(bottom=0.2)

    
    plt.plot(indices,angles,label="Angle between B and K")
    delta = int(max(angles) - min(angles))
    offs = min(angles)
    plt.plot((index_reversal)*np.ones(delta),offs+np.array(range(delta)),'y',label="Time/Location of Faraday Rotation Inflection")
    plt.plot((indx_closest)*np.ones(delta),offs+np.array(range(delta)),'g',label="Approximate Location of closest approach (" + str(dists[indx_closest]) + " km)")
    plt.title("Relative angle of B vector vs. K vector for CASSIOPE ephemeris on " + str(date_string))
    #plt.xlabel('Time elapsed during pass (seconds)')
    plt.xticks(tick_indices, my_xticks)
    plt.legend()
    plt.show()

def plot_kvec(date_string):
    """
    A function for visualizing the k-vector's time evolution throughout a 
    pass of EPOP.

    *** PARAMS ***
    date_string (string): string in form "20160418" to denote date for which to plot KB angle

    *** RETURNS ***
    - (just plots)

    """
    fname,index_reversal = get_ottawa_data(date_string)
    lons,lats,alts,ephtimes,mlons,mlats,mlts,pitch,yaw,roll = get_rri_ephemeris_full(fname)
    bvecs,kvecs,angles = get_kb_ottawa_angle(lons,lats,alts,ephtimes)
    indx_closest, dists = get_closest_ottawa_approach(lons,lats,alts)
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
    dist_t = dists[0]
    my_xticks.append("Time (UTC):    "+str(dt_t.time())+"\nLatitude:    "+str(lat_t)+\
        "\nLongitude:    "+str(lon_t)+"\nAltitude:    "+str(alt_t)+\
    "\nMagnetic Local Time:    "+str(mlt_t)+"\nMagnetic Latitude:    "+str(mlat_t)+\
    "\nMagnetic Longitude:    "+str(mlon_t)+"\nDistance (km):    "+str(dist_t))
    for i in range(num_ticks-1):
        alt_t = alts[tick_sep*(i+1)]
        lon_t = lons[tick_sep*(i+1)]
        lat_t = lats[tick_sep*(i+1)]
        mlon_t = mlons[tick_sep*(i+1)]
        mlat_t = mlats[tick_sep*(i+1)]
        dt_t  = times[tick_sep*(i+1)]
        mlt_t = mlts[tick_sep*(i+1)]
        dist_t = dists[tick_sep*(i+1)]
        my_xticks.append(str(dt_t.time())+"\n"+str(lat_t)+"\n"+str(lon_t)+"\n"+str(alt_t)+\
                "\n"+str(mlt_t)+"\n"+str(mlat_t)+"\n"+str(mlon_t)+"\n"+str(dist_t))
    
    indices = range(angles.__len__())
    tick_indices = [i*tick_sep for i in range(num_ticks)]
    
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = Axes3D(fig)
    kx = [kv[0] for kv in kvecs]
    ky = [kv[1] for kv in kvecs]
    kz = [kv[2] for kv in kvecs]
    Axes3D.plot(ax,kx,ky,zs=kz)
    plt.xlabel("Kx direction (North)")
    plt.ylabel("Ky direction (East)")
    ax.set_zlabel("Kz direction (Down)")
    plt.title("Change of K vector components during RRI pass on " + str(date_string))
    plt.show() 
 
def plot_ramdir(date_string):
    """
    Makes a 3D plot of the ram direction components at each ephemeris point for
    the given date.  
    """
    lons,lats,alts,ephtimes,mlons,mlats,mlts,pitch,yaw,roll = get_rri_ephemeris_full(datname)
    vs,dists = get_ramdirs(lons,lats,alts,ephtimes)   
    indx_closest, dists = get_closest_ottawa_approach(lons,lats,alts)
    times = ephems_to_datetime(ephtimes)
    
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = Axes3D(fig)
    vx = [v[0] for v in vs]
    vy = [v[1] for v in vs]
    vz = [v[2] for v in vs]
    Axes3D.plot(ax,vx,vy,zs=vz)
    plt.xlabel("Vx direction - North (km/s)")
    plt.ylabel("Vy direction - East (km/s)")
    ax.set_zlabel("-Vz direction - Up (km/s)")
    plt.title("Change of V vector components during RRI pass on " + str(date_string))
    plt.show() 

def plot_kdip_angle(date_string):
    """
    Plots the angle between the k_LOS vector and the dip_dir vector. 
    If they were perfectly matched, we would expect exactly 180 degrees.
    """
    fname,index_reversal = get_ottawa_data(date_string)
    lons,lats,alts,ephtimes,mlons,mlats,mlts,pitch,yaw,roll = get_rri_ephemeris_full(fname)
    dipole_dirs, kdip_angles = get_kdip_angles(lons,lats,alts,ephtimes,pitch,yaw,roll)
    #plot_kb_angle(date_string)
    #plot_kvec(date_string)
    #plot_ramdir("20160418")
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.subplots_adjust(bottom=0.2)
    plt.plot(kdip_angles)
    ephem_ticks(lons,lats,alts,ephtimes,mlons,mlats,mlts)
    plt.ylabel('Angle (degrees)')
    plt.title('Plot of angle between K_los from Ottawa transmitter and Dipole direction of RRI for ' + date_string)
    plt.show()
    

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]]) 

def get_kdip_angles(lons,lats,alts,ephtimes,pitch,yaw,roll):
    """
    Call this function to retrieve a list of the direction of the RRI dipole plane 
    in N-E-down coordinates for the first n-1 ephemeris points provided, as well
    as the relative angle between the K_los vector from the Ottawa transmitter
    and the direction of the RRI dipole plane.

    """
    # vs is a vector of ram directions in N, E, Down coordinates
    vs,dists = get_ramdirs(lons,lats,alts,ephtimes)
   
    # kvecs will be k_LOS vectors in N, E, Down coordinates 
    kvecs = []
    for i in range(lons.__len__()):
        kvec = get_ottawa_kvec(lons[i],lats[i],alts[i],ephtimes[i])
        kvecs.append(kvec)

    # Word of Gareth:
    # x is ram direction, z is nadir direction, y is Z cross X.
    xdirs = vs
    zdirs = np.array([ (0,0,1) for i in range(xdirs.__len__())])
    ydirs = np.cross(zdirs,xdirs)
    
    # yaw: rot around z, pitch: rot around y, roll: rot around x
    # Assuming as seems to be confirmed in documentation that the Dipole is in the
    # x-direction on CASSIOPE, and thus, in default position, towards the ram direction.
    dipole_dirs = []
    kdip_angles = []
    for i in range(xdirs.__len__()):
        yaw_rot = rotation_matrix(zdirs[i],np.deg2rad(yaw[i]))
        pitch_rot = rotation_matrix(ydirs[i],np.deg2rad(pitch[i]))
        roll_rot = rotation_matrix(xdirs[i],np.deg2rad(roll[i]))
        initial_dipole_vec = xdirs[i]
        intermed1 = np.dot(yaw_rot,initial_dipole_vec)
        intermed2 = np.dot(pitch_rot,intermed1)
        # After all the rotations, we have dip_dir
        dip_dir = np.dot(roll_rot,intermed2)
        dipole_dirs.append(dip_dir)
        # Determine what the relative angle is between k_LOS and dip_dir
        kdip_angle = np.arccos(np.dot(dip_dir,kvecs[i])/(np.linalg.norm(dip_dir)*np.linalg.norm(kvecs[i])))
        kdip_angles.append(np.rad2deg(kdip_angle))
    return dipole_dirs, kdip_angles    

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
   

def plot_1822():
    """
    Plot April 18th and 22nd ephemeris tracks on same mapobj
    """
    # TODO: fixup this documentation
    date_string="20160418"
       
    fname,index_reversal,ellip_reversal = get_ottawa_data2(date_string)
    geog_longs,geog_lats,alts,ephemtimes = get_rri_ephemeris(fname)

    fname_22nd, idx_rev_22nd, ellip_reversal_22nd = get_ottawa_data2("20160422")
    lons_22nd, lats_22nd, alts_22nd, ephtimes_22nd = get_rri_ephemeris(fname_22nd)

    times = ephems_to_datetime(ephemtimes)
    times_22nd = ephems_to_datetime(ephtimes_22nd)

    indx_shortest, dists = get_closest_ottawa_approach(geog_longs, geog_lats, alts)
    appr_time = times[indx_shortest]

    indx_shortest_22nd, dists_22nd = get_closest_ottawa_approach(lons_22nd, lats_22nd, alts_22nd)
    appr_time_22nd = times_22nd[indx_shortest_22nd]
    
    # The numeric data type that I was retrieving from geog_longs, when _NOT_ stored
    # in an array, was being rejected by the mapObj() function below. So I convert 
    # these numbers to floats explicitly here.
    shortest_ephem_long = float(geog_longs[indx_shortest])
    shortest_ephem_lat = float(geog_lats[indx_shortest])
    inversion_ephem_long = float(geog_longs[index_reversal])
    inversion_ephem_lat = float(geog_lats[index_reversal])
    elliplon_18th = float(geog_longs[ellip_reversal])
    elliplat_18th = float(geog_lats[ellip_reversal])

    elliplon_22nd = float(lons_22nd[ellip_reversal_22nd])
    elliplat_22nd = float(lats_22nd[ellip_reversal_22nd])
    shlon_22nd = float(lons_22nd[indx_shortest])
    shlat_22nd = float(lats_22nd[indx_shortest])
    invlon_22nd = float(lons_22nd[index_reversal])
    invlat_22nd = float(lats_22nd[index_reversal])
    
    # A different font for the legend etc. might be nice
    fig = plt.figure(1, figsize=(14,8))
    font = {'fontname':'Computer Modern'}
    #m = plotUtils.mapObj(lat_0=45.0, lon_0=-75.0, width=60e3*180, height=50e3*90, coords='geo', resolution='i',datetime=times[0])
    m = plotUtils.mapObj(lat_0=57.0, lon_0=0.0, width=25e3*180, height=25e3*90, coords='mag', resolution='i', datetime=times[0])

    # FIRST: Plot the location of Ottawa
    x,y = m(OTTAWA_TX_LON,OTTAWA_TX_LAT,coords='geo')
    m.plot(x,y,'r-o',markersize=9,label="Ottawa")
    
    # SECOND: Plot the satellite ground-track.
    # Day of 18th
    x,y = m(geog_longs, geog_lats, coords='geo')
    m.plot(x,y,'b',label="EPOP ground track")
    # Day of 22nd
    x,y = m(lons_22nd, lats_22nd, coords='geo')
    m.plot(x,y,'b')
    
    # THIRD: Plot a circle emphasizing the point of closest approach
    x,y = m(shortest_ephem_long,shortest_ephem_lat, coords='geo')
    m.plot(x,y,'bo')
    x,y = m(shlon_22nd, shlat_22nd, coords='geo')
    m.plot(x,y,'bo',label=("Closest Approach TX"))

    # FOURTH: Plot the line from Ottawa to the nearest approach of the satellite.
    x,y = m([shortest_ephem_long, OTTAWA_TX_LON], [shortest_ephem_lat, OTTAWA_TX_LAT], coords='geo')
    m.plot(x,y,'g')
    x,y = m([shlon_22nd, OTTAWA_TX_LON], [shlat_22nd, OTTAWA_TX_LAT], coords='geo')
    m.plot(x,y,'g')
    
    
    # FIFTH: Plot the piont I've determined is the point of the Faraday Rotation inversion.
    x,y = m(inversion_ephem_long,inversion_ephem_lat,coords='geo')
    m.plot(x,y,'yo',label=("Inflection of Faraday Rotation"))
    x,y = m(invlon_22nd, invlat_22nd, coords='geo')
    m.plot(x,y,'yo')

    # SIXTH: Plot the point we've determined to be the ellipticity reversal
    x,y = m(elliplon_18th, elliplat_18th, coords='geo')
    m.plot(x,y,'go',markersize=4,label=("Inflection of Ellipticity Angle"))
    x,y = m(elliplon_22nd, elliplat_22nd, coords='geo')
    m.plot(x,y,'go',markersize=4)

    # SEVENTH: GET IGRF DATA FOR EACH EPHEMERIS POINT
    """
    itype = 1 #Geodetic coordinates
    pyDate = times[0] # The first time we pull from the RRI file
    date = utils.dateToDecYear(pyDate) # decimal year
    alt = 300. # altitude #TODO: grab altitudes for series of satellite positions we care about.
    stp = 1. #
    xlti, xltf, xltd = OTTAWA_TX_LAT, OTTAWA_TX_LAT,stp # latitude start, stop, step
    xlni, xlnf, xlnd = OTTAWA_TX_LON, OTTAWA_TX_LON,stp # longitude start, stop, step
    ifl = 0 # Main field
    # Call fortran subroutine
    lat,lon,d,s,h,x,y,z,f = igrf.igrf11(itype,date,alt,ifl,xlti,xltf,xltd,xlni,xlnf,xlnd)
    """
    ax = plt.gca()
    ax.set_xlabel('Magnetic Longitude (degrees)')
    ax.xaxis.set_label_coords(0.5,-0.050)
    #plt.xlabel('Magnetic Longitude (degrees)')
    plt.ylabel('Magnetic Latitude (degrees)')
    plt.title("ePOP Pass vs. Ottawa radar on April 18th, 2016 (left) and April 22nd, 2016 (right)")
    plt.legend(loc='lower right', numpoints = 1)#loc='best')
    
    plt.savefig('tmp_1822.eps', format='eps', bbox_inches='tight')
    plt.savefig('tmp_1822.png', format='png', bbox_inches='tight')
    plt.show()

def plot_all5():
    """
    Plot April 18th-22nd ephemeris tracks on same mapobj
    """
    fname_18th, idx_rev_18th, ellip_rev_18th = get_ottawa_data2("20160418")
    lons_18th, lats_18th, alts_18th, ephtimes_18th = get_rri_ephemeris(fname_18th)
    fname_19th, idx_rev_19th, ellip_rev_19th = get_ottawa_data2("20160419")
    lons_19th, lats_19th, alts_19th, ephtimes_19th = get_rri_ephemeris(fname_19th)
    fname_20th, idx_rev_20th, ellip_rev_20th = get_ottawa_data2("20160420")
    lons_20th, lats_20th, alts_20th, ephtimes_20th = get_rri_ephemeris(fname_20th)
    fname_21st, idx_rev_21st, ellip_rev_21st = get_ottawa_data2("20160421")
    lons_21st, lats_21st, alts_21st, ephtimes_21st = get_rri_ephemeris(fname_21st)
    fname_22nd, idx_rev_22nd, ellip_rev_22nd = get_ottawa_data2("20160422")
    lons_22nd, lats_22nd, alts_22nd, ephtimes_22nd = get_rri_ephemeris(fname_22nd)
    
    times_18th = ephems_to_datetime(ephtimes_18th)
    times_19th = ephems_to_datetime(ephtimes_19th)
    times_20th = ephems_to_datetime(ephtimes_20th)
    times_21st = ephems_to_datetime(ephtimes_21st)
    times_22nd = ephems_to_datetime(ephtimes_22nd)

    indx_shortest_18th, dists_18th = get_closest_ottawa_approach(lons_18th, lats_18th, alts_18th)
    indx_shortest_19th, dists_19th = get_closest_ottawa_approach(lons_19th, lats_19th, alts_19th)
    indx_shortest_20th, dists_20th = get_closest_ottawa_approach(lons_20th, lats_20th, alts_20th)
    indx_shortest_21st, dists_21st = get_closest_ottawa_approach(lons_21st, lats_21st, alts_21st)
    indx_shortest_22nd, dists_22nd = get_closest_ottawa_approach(lons_22nd, lats_22nd, alts_22nd)
    
    # The numeric data type that I was retrieving from geog_longs, when _NOT_ stored
    # in an array, was being rejected by the mapObj() function below. So I convert 
    # these numbers to floats explicitly here.
    shlon_18th = float(lons_18th[indx_shortest_18th])
    shlat_18th = float(lats_18th[indx_shortest_18th])
    invlon_18th = float(lons_18th[idx_rev_18th])
    invlat_18th = float(lats_18th[idx_rev_18th])
    elliplon_18th = float(lons_18th[ellip_rev_18th])
    elliplat_18th = float(lats_18th[ellip_rev_18th])
 
    shlon_19th = float(lons_19th[indx_shortest_19th])
    shlat_19th = float(lats_19th[indx_shortest_19th])
    invlon_19th = float(lons_19th[idx_rev_19th])
    invlat_19th = float(lats_19th[idx_rev_19th])
    elliplon_19th = float(lons_19th[ellip_rev_19th])
    elliplat_19th = float(lats_19th[ellip_rev_19th])

    shlon_20th = float(lons_20th[indx_shortest_20th])
    shlat_20th = float(lats_20th[indx_shortest_20th])
    invlon_20th = float(lons_20th[idx_rev_20th])
    invlat_20th = float(lats_20th[idx_rev_20th])
    elliplon_20th = float(lons_20th[ellip_rev_20th])
    elliplat_20th = float(lats_20th[ellip_rev_20th])

    shlon_21st = float(lons_21st[indx_shortest_21st])
    shlat_21st = float(lats_21st[indx_shortest_21st])
    invlon_21st = float(lons_21st[idx_rev_21st])
    invlat_21st = float(lats_21st[idx_rev_21st])
    elliplon_21st = float(lons_21st[ellip_rev_21st])
    elliplat_21st = float(lats_21st[ellip_rev_21st])

    shlon_22nd = float(lons_22nd[indx_shortest_22nd])
    shlat_22nd = float(lats_22nd[indx_shortest_22nd])
    invlon_22nd = float(lons_22nd[idx_rev_22nd])
    invlat_22nd = float(lats_22nd[idx_rev_22nd])
    elliplon_22nd = float(lons_22nd[ellip_rev_22nd])
    elliplat_22nd = float(lats_22nd[ellip_rev_22nd])
   
    # A different font for the legend etc. might be nice
    fig = plt.figure(1, figsize=(14,8))
    font = {'fontname':'Computer Modern'}
    #m = plotUtils.mapObj(lat_0=45.0, lon_0=-75.0, width=50e3*180, height=40e3*90, coords='geo',resolution='i',datetime=times_20th[0])
    m = plotUtils.mapObj(lat_0=57.0, lon_0=0.0, width=25e3*180, height=25e3*90, coords='mag',resolution='i',datetime=times_20th[0])

    # FIRST: Plot the location of Ottawa
    x,y = m(OTTAWA_TX_LON,OTTAWA_TX_LAT,coords='geo')
    m.plot(x,y,'r-o',markersize=9,label="Ottawa")
    
    # SECOND: Plot the satellite ground-track.
    x,y = m(lons_18th, lats_18th, coords='geo')
    m.plot(x,y,'b-',label="April 18th")#label="EPOP ground track")
    x,y = m(lons_19th, lats_19th, coords='geo')
    m.plot(x,y,'b:',label="April 19th")
    x,y = m(lons_20th, lats_20th, coords='geo')
    m.plot(x,y,'b--',label="April 20th")
    x,y = m(lons_21st, lats_21st, coords='geo')
    m.plot(x,y,'k-',label="April 21st")
    x,y = m(lons_22nd, lats_22nd, coords='geo')
    m.plot(x,y,'k:',label="April 22nd")
    
    # THIRD: Plot a circle emphasizing the point of closest approach
    x,y = m(shlon_18th, shlat_18th, coords='geo')
    m.plot(x,y,'bo',label='Closest Approach TX')
    x,y = m(shlon_19th, shlat_19th, coords='geo')
    m.plot(x,y,'bo')
    x,y = m(shlon_20th, shlat_20th, coords='geo')
    m.plot(x,y,'bo')
    x,y = m(shlon_21st, shlat_21st, coords='geo')
    m.plot(x,y,'bo')
    x,y = m(shlon_22nd, shlat_22nd, coords='geo')
    m.plot(x,y,'bo')

    # FOURTH: Plot the piont I've determined is the point of the Faraday Rotation inversion.
    x,y = m(invlon_18th, invlat_18th, coords='geo')
    m.plot(x,y,'yo',markersize=6,label=("Inflection of Faraday Rotation"))
    x,y = m(invlon_19th, invlat_19th, coords='geo')
    m.plot(x,y,'yo')
    x,y = m(invlon_20th, invlat_20th, coords='geo')
    m.plot(x,y,'yo')
    x,y = m(invlon_21st, invlat_21st, coords='geo')
    m.plot(x,y,'yo')
    x,y = m(invlon_22nd, invlat_22nd, coords='geo')
    m.plot(x,y,'yo')

    # FIFTH: Plot the point we've determined to be the ellipticity reversal
    x,y = m(elliplon_18th, elliplat_18th, coords='geo')
    m.plot(x,y,'go',markersize=4,label=("Inflection of Ellipticity Angle"))
    x,y = m(elliplon_19th, elliplat_19th, coords='geo')
    m.plot(x,y,'go',markersize=4)
    x,y = m(elliplon_20th, elliplat_20th, coords='geo')
    m.plot(x,y,'go',markersize=4)
    x,y = m(elliplon_21st, elliplat_21st, coords='geo')
    m.plot(x,y,'go',markersize=4)
    x,y = m(elliplon_22nd, elliplat_22nd, coords='geo')
    m.plot(x,y,'go',markersize=4)

    # FIFTH: a few lines of magnetic longitude and latitude will be plotted as well.
    #N = 10
    #latdivs = 90./N
    #londivs = 360./N
    #for n in range(N):
    #    merid_mlat = np.arange(181) - 90.
    #    merid_mlon = merid_mlat*0 - n*londivs
    #    paral_mlon = np.arange(361) - 180.
    #    paral_mlat = paral_mlon*0. + n*latdivs
    #    x,y = m(merid_mlon, merid_mlat, coords='mag')
    #    if n == 0:
    #        m.plot(x,y,'k',label="Magnetic Longitude and Latitude")
    #        continue
    #    m.plot(x,y,'k')#,label="Line of Magnetic Longitude of -20 Degrees")
    #    x,y = m(paral_mlon, paral_mlat, coords='mag')
    #    m.plot(x,y,'k')#,label="Line of Magnetic Latitude of +75 Degrees")
    ax = plt.gca()
    ax.set_xlabel('Magnetic Longitude (degrees)')
    ax.xaxis.set_label_coords(0.5,-0.050)
    #plt.xlabel('Magnetic Longitude (degrees)')
    plt.ylabel('Magnetic Latitude (degrees)')
    plt.title("ePOP Pass vs. Ottawa radar from April 18th - 22nd, 2016")
    plt.legend(loc='lower right', numpoints = 1)#loc='best')
    plt.savefig('tmp_all5.eps', format='eps', bbox_inches='tight')
    plt.savefig('tmp_all5.png', format='png', bbox_inches='tight')
    plt.show()
 
 

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    lon = OTTAWA_TX_LON
    lat = OTTAWA_TX_LAT
    alt = OTTAWA_TX_ELEV 
    #(bx,by,bz) = get_bvec(lon,lat,alt,dt.now())
    
    date_string = "20160422"
    datpath,datname = initialize_data()
    
    fname,index_reversal = get_ottawa_data(date_string)
    lons,lats,alts,ephtimes,mlons,mlats,mlts,pitch,yaw,roll = get_rri_ephemeris_full(fname)
    dipole_dirs, kdip_angles = get_kdip_angles(lons,lats,alts,ephtimes,pitch,yaw,roll)
    
    #plot_kb_angle(date_string)
    #plot_kvec(date_string)
    #plot_ramdir("20160418")
    #plot_kdip_angle("20160418")
    """
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.subplots_adjust(bottom=0.2)
    plt.plot(kdip_angles)
    ephem_ticks(lons,lats,alts,ephtimes,mlons,mlats,mlts)
    plt.ylabel('Angle (degrees)')
    plt.title('Plot of angle between K_los from Ottawa transmitter and Dipole direction of RRI for ' + date_string)
    plt.show()
    """
    
