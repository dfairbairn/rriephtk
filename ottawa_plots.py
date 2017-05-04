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
    Takes two lon/lat/alt points and determines the elevation angle necessary
    for getting from the first point to the second.

    ***PARAMS***
        lon1 [float]: longitude of point 1 (deg)
        lat1 [float]: latitude of point 1 (deg)
        alt1 [float]: altitude of opint 1 (km)
        lon2 [float]: longitude of point 2 (deg)
        lat2 [float]: latitude of point 2 (deg)
        alt2 [float]: altitude of point 2 (km)

    ***RETURNS***
        elev_angle [float]: elevation angle 

    *currently doesn't account for Earth's curvature. 
 
    """
    arcdist = haversine(lon1, lat1, lon2, lat2)
    delta_alt = alt2 - alt1
    # doesn't account for curvature of earth right now! #TODO: curvature 
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
    
