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

OTTAWA_TX_LON = -75.552
OTTAWA_TX_LAT = 45.403
def get_ottawa_kvec(glon, glat, altitude, time):
    """
    This function takes a satellite ephemeris point as input, calculating
    the straight-line k-vector from the Ottawa-based transmitter in terms of 
    North (x), East (y), Down (z) components.
    """
    # The specific coordinates of the NRCAN geomagnetic    
    OTTAWA_TX_LON = -75.552 
    OTTAWA_TX_LAT = 45.403

    # In spherical coordinates, subtracting vectors doesn't get us the path
    # from one point to another along the surface. To accomplish that, we use
    # bearings and elevation angles.
    
    init_bearing,final_bearing = get_bearing(OTTAWA_TX_LON, OTTAWA_TX_LAT, glon, glat)
    elev_angle = get_elevation_angle(OTTAWA_TX_LON, OTTAWA_TX_LAT, 0., glon, glat, altitude)

    # for calculation of vector components, we need theta (the bearing in the
    # N-E plane) at CASSIOPE, and phi, the angle-from-down (toward the N-E plane)
    theta = final_bearing
    phi = 90. + elev_angle
    kx = np.sin(np.deg2rad(phi))*np.cos(np.deg2rad(theta))
    ky = np.sin(np.deg2rad(phi))*np.sin(np.deg2rad(theta))
    kz = np.cos(np.deg2rad(theta)) 
    return np.array((kx,ky,kz))

def get_bearing(lon1, lat1, lon2, lat2):
    """
    Gives the bearing clockwise from north in degrees from point 1 to point2,
    at point 1 (init_bearing) and at point 2 (final_bearing)
    """
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

def get_nvec(lon, lat, alt):
    """
    This function calculates and returns an n-vector as described at the 
    website: http://www.movable-type.co.uk/scripts/latlong-vectors.html

    **PARAMS*
    lon (float): 
    lat (float):
    alt (float): altitude of the point

    *** MAY NOT BE USEFUL FOR THINGS IM DOING AS OF AUG 11 ***
    """
    latrad = np.deg2rad(lat)
    lonrad = np.deg2rad(lon)
    vx = np.cos(latrad)*np.cos(lonrad)
    vy = np.cos(latrad)*np.sin(lonrad)
    vz = np.sin(latrad)
    return vx,vy,vz

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

def get_closest_ottawa_approach(glons, glats):
    """
    Uses a haversine formula-based distance calculation to find the closest 
    approach of a set of points describing a ground-track.

    *** PARAMS *** 
    glons (np.array of floats): longitudes of ground-track points
    glats (np.array of floats): latitudes of ground-track points

    ** RETURNS **
    indx_shortest (integer): index of the lat/lon pair where the ottawa transmitter is closest
    dist_shortest (float): the actual closest distance
    """
    #TODO: Note: the 'closest distance' doesn't account for altitude! Do this?

    # *** FINDING THE CLOSEST APPROACH ***
    # Using the Haversine formula (in a function in script_utils.py), the closest
    # approach is determined by brute force.
    dists = []
    longdists = []
    latdists = []
    dist_shortest = sys.maxint #Initially set this very high
    for i in range(np.size(glons)):
        dist = haversine(glons[i], glats[i], OTTAWA_TX_LON, OTTAWA_TX_LAT)
        # Initially I took a quick and dirty approach to find the point of smallest
        # Euclidean distance in terms of latitudes and longitudes.
        #longdist = abs(geog_longs[i] - OTTAWA_TX_LON) # difference of longitudes
        #latdist = abs(geog_lats[i] - OTTAWA_TX_LAT) # difference of latitudes
        #dist = np.sqrt(longdist*longdist + latdist*latdist)  
        if dist < dist_shortest:
            dist_shortest = dist
            indx_shortest = i
    return indx_shortest, dist_shortest

def get_ottawa_data(date_string):
    """
    Right now, this just returns the ephemeris data from CASSIOPE on this date,
    as well as a hardcoded-by-inspection index/second number that I determined, 
    but could be extended to do analytically in the future. 
   
    *** PARAMS ***
    date_string (string): a string of the form "20160418" to denote april 18th of 2016 

    *** RETURNS ***
    glons
    glats
    alts
    ephemtimes
    index_reversal (integer): 
    """
    if isinstance(date_string, type(None)): date_string="20160418"

    # **** CHOOSE ONE OF THESE RRI FILES THEN RUN THE SCRIPT ****
    if "20160418"==date_string:
        geog_longs,geog_lats,alts,ephemtimes = get_rri_ephemeris("./data/RRI_20160418_222759_223156_lv1_v2.h5") #18th
        index_reversal = 167 #for 18th
    elif "20160419"==date_string:
        index_reversal = 178 #for 19th
        geog_longs,geog_lats,alts,ephemtimes = get_rri_ephemeris("./data/RRI_20160419_220939_221336_lv1_v2.h5") #19th
    elif "20160420"==date_string:
        index_reversal = 213 #for 20th
        geog_longs,geog_lats,alts,ephemtimes = get_rri_ephemeris("./data/RRI_20160420_215117_215514_lv1_v2.h5") #20th
    elif "20160421"==date_string:
        index_reversal = 205 #?? for 21st?
        geog_longs,geog_lats,alts,ephemtimes = get_rri_ephemeris("./data/RRI_20160421_213255_213652_lv1_v2.h5") #21st
    elif "20160422"==date_string:
        index_reversal = 222 #for 22nd
        geog_longs,geog_lats,alts,ephemtimes = get_rri_ephemeris("./data/RRI_20160422_211435_211832_lv1_v2.h5") #22nd
    else:
        print "Invalid input date."
        return None    
    return geog_longs,geog_lats,alts,ephemtimes,index_reversal
   


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

    """
    # TODO: fixup this documentation
    if isinstance(date_string, type(None)): date_string="20160418"
       
    geog_longs,geog_lats,alts,ephemtimes,index_reversal = get_ottawa_data(date_string)

    # Location of Ottawa: I looked it up and hard-coded it at the top
    times = ephems_to_datetime(ephemtimes)

    indx_shortest, dist_shortest = get_closest_ottawa_approach(geog_longs, geog_lats)
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
    m.plot(x,y,'ro',label="Ottawa")
    
    # SECOND: Plot the satellite ground-track.
    x,y = m(geog_longs, geog_lats, coords='geo')
    m.plot(x,y,'b',label="EPOP ground track")
    
    # THIRD: Plot a circle emphasizing the point of closest approach
    x,y = m(shortest_ephem_long,shortest_ephem_lat, coords='geo')
    m.plot(x,y,'bo',label=("Shortest Approach at " + str(appr_time)))
    
    # FOURTH: Plot the line from Ottawa to the nearest approach of the satellite.
    x,y = m([shortest_ephem_long, OTTAWA_TX_LON], [shortest_ephem_lat, OTTAWA_TX_LAT], coords='geo')
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
    xlti, xltf, xltd = OTTAWA_TX_LAT, OTTAWA_TX_LAT,stp # latitude start, stop, step
    xlni, xlnf, xlnd = OTTAWA_TX_LON, OTTAWA_TX_LON,stp # longitude start, stop, step
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
OTTAWA_TX_LON = -75.552
OTTAWA_TX_LAT = 45.403
lon = OTTAWA_TX_LON
lat = OTTAWA_TX_LAT
alt = 0.
#(bx,by,bz) = get_bvec(lon,lat,alt,dt.now())

date_string = "20160418"
datpath,datname = initialize_data()
lons,lats,alts,times,index_reversal = get_ottawa_data(date_string)
bvecs,kvecs,angles = get_kb_ottawa_angle(lons,lats,alts,times)
indx_closest, dist_closest = get_closest_ottawa_approach(lons,lats)
plt.plot(angles,label="Angle between B and K")
plt.plot((index_reversal)*np.ones(100),range(100),'y',label="Time/Location of Faraday Rotation Reversal")
plt.plot((indx_closest)*np.ones(100),range(100),'g',label="Approximate Location of closest approach")
plt.title("Relative angle of B vector vs. K vector for CASSIOPE ephemeris on " + str(date_string))
plt.xlabel('Time elapsed during pass (seconds)')
plt.legend()
plt.show()
