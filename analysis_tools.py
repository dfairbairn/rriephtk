# -*- coding: utf-8 -*-
"""
file: 'analysis_tools.py'
description:
    This file bundles together the variety of tools developed for calculating 
    factors associated with RRI experiments. 

    These include telemetry-related functions, mathematical transformation
    functions, parameter calculation, etc.

author: David Fairbairn
date: May 2017

"""
import sys
import numpy as np
import matplotlib.pyplot as plt

import data_utils
import magnet_data

from script_utils import *

OTTAWA_TX_LON = -75.552
OTTAWA_TX_LAT = 45.403
OTTAWA_TX_ELEV = 0.070 # 70m elevation in Ottawa - small compared to 300km satellite

# ------------------------------------------------------------------------------
#                       Index of Refraction-related Code
#                       --------------------------------
#    The purpose of this code is the calculation of index of refraction vs
#    incidence angle for radio wave propagation direction compared to the
#    ionospheric plasma's orientation/local B field
# ------------------------------------------------------------------------------

def basic_ionosphere_params(altitude=300.):
    """ Returns cyclotron frequency, debye length, and plasma frequency in 
    ionosphere for a given altitude. 
  
    *PARAMS*:
        [altitude]: In km. Must be in ionosphere for this model to work. default 300.
 
    ANGULAR FREQUENCIES BTW

    *RETURNS*:
        omega_p: Plasma frequency (angular)
        omega_c: Cyclotron frequency (angular)
        l_d: Debye length
     
    Based on:
    https://smallsats.org/2013/03/25/debye-length-plasma-frequency-and-cyclotron-frequency/
    """

    if altitude > 400. or altitude < 200.:
        print "Bad altitude parameter"
        return None

    # Assumed parameters:
    me   = 9.109E-31            #[kg]           Electron rest mass
    mp   = 1.673E-27            #[kg]           Proton rest mass
    eps0 = 8.8542E-12           #[A*s/(V*m)]    Permittivity
    e    = 1.602E-19            #[C]            Elementary charge
    Re   = 6.37E6               #[m]            Earth’s radius
    Md   = 8E15                 #[T*m3]         Earth’s magnetic dipole moment

    n   = 1E12                     #[m-3]  Electron density < HIGH >
    kTe = 0.1                      #[eV]   Thermal energy
    r   = Re + altitude*1000.      #[km]   Radial distance from the Earth’s center
    omega_p  = ((e**2)*n/(eps0*me))**0.5  #[Hz]   Plasma frequency
    l_d  = (eps0*kTe/(n*e))**0.5     #[m]    Debye length
    B   = 2*Md/(r**3)                #[T]    Earth’s magnetic field
    omega_c  = e*B/(me)            #[Hz]   Cyclotron frequency

    return omega_c,omega_p,l_d

def plasma_freq(n_e):
    """
    Given an electron density parameter (n_e), compute the plasma frequency.
    """
    eps0 = 8.8542E-12           #[A*s/(V*m)]    Permittivity
    e    = 1.602E-19            #[C]            Elementary charge
    me   = 9.109E-31            #[kg]           Electron rest mass
    omega_p  = ((e**2)*n_e/(eps0*me))**0.5  #[Hz]   Plasma frequency
    return omega_p

def cyclotron_freq(B):
    """
    Takes a magnetic field value(s) parameter and uses it to calculate the
    cyclotron frequency for an electron in this field.

    Basically either uses a single B field for a cyclotron frequency, or takes the 
    average of an array and uses that to find a cyclotron frequency.

    ***PARAMS***
        B [float] or [list/ndarray of floats]: magnetic field. If in list form,
            this function takes the average magnitude of B's [nT]
         *Note that B is expected to be in nanotesla (e.g. around 50000nT)

    ***RETURNS***
        omega_c [float]: cyclotron frequency (rad/s)
    """
    e    = 1.602E-19            #[C]            Elementary charge
    me   = 9.109E-31            #[kg]           Electron rest mass
    if type(B)==np.ndarray or type(B)==list:
        B_mag = 1E-9*np.mean( [np.linalg.norm(j) for j in B ]) # 1E-9 for m vs km ^3
    elif isinstance(B, float):
        B_mag = B
    else: 
        print "Unexpected data type for B field parameter"
        print (B,"\t",type(B))
        return None
    omega_c  = e*B_mag/(me)            #[Hz]   Cyclotron frequency
    return omega_c

def improved_cyclotron_freq(lons,lats,alts,ephtimes,fof2_alt=250.):
    """

    """
    import magnet_data 
    e =  1.602E-19            #[C]            Elementary charge
    me = 9.109E-31            #[kg]           Electron rest mass

    # attempt to find B field within peak plasma density region at fof2 altitude
    B = magnet_data.get_igrf(lons,lats,fof2_alt*np.ones(len(alts)),ephtimes) 
    B_mag = 1E-9*np.mean( [np.linalg.norm(j) for j in B ]) # 1E-9 for m vs km ^3
    omega_c = e*B_mag/me
    return omega_c
    

def appleton_coeffs(omega_c, omega_p, freq, nu_e):
    """
    Gets the Appleton Hartree X,Y,Z coefficients given the ionospheric params.

    Parameters:
        omega_c: angular cyclotron frequency (rad/s)
        omega_p: angular plasma frequency (rad/s)
        omega:  angular frequency of radio wave (rad/s)
        nu_e: electron collision frequency (Hz)

    Returns:
        X = (omega_p/omega)^2
        Y = (omega_c/omega)
        Z = nu/omega

    """
    omega = (2*np.pi)*freq
    X = (omega_p/omega)**2
    Y = omega_c/omega
    Z = nu_e/omega
    return X,Y,Z


def appleton_hartree(X,Y,Z,theta):
    """
    Computes the index of refraction of the ionospheric plasma medium
    using the appleton-hartree equation. Aided by helper functions for X, Y, 
    and Z.
    
    Parameters: 
        X: X parameter in Appleton Hartree
        Y: Y parameter in Appleton Hartree
        Z: Z parameter in Appleton Hartree
        theta_i: aspect angles (potentially multiple) in radians

    Returns:
        nplus: The positive index of refraction from Appleton-Hartree
        nminus: The negative index of refraction from Appleton-Hartree
    """

    term1 = 1. - 1j*Z
    term2 = ( (Y*np.sin(theta))**2)/(2*(1 - X - 1j*Z))
    term3 = np.sqrt( ( ( Y*np.sin(theta) )**4 )/(4*( (1 - X - 1j*Z)**2 )) + (Y*np.cos(theta))**2 )
    nplus_squared = 1 - X/(term1 - term2 + term3)
    nminus_squared = 1 - X/(term1 - term2 - term3)
    nplus = np.sqrt(nplus_squared)
    nminus = np.sqrt(nminus_squared)
    return nplus,nminus

def txpass_to_indices(lons,lats,alts,ephtimes,tx_lon,tx_lat,freq,omega_p,omega_c):
    """
    Calculates the expected indices of refraction at the fof2 peak for radio
    waves directed from a specified transmitter toward CASSIOPE.

    Each satellite ephemeris point ([lon,lat,alt,ephtime]) receives aspect
    angles (from get_kb_angle) and corresponding indices of refraction for the
    ionospheric plasma at the given frequency.

    ***PARAMS***
        lons [float]: longitude of satellite (deg)
        lats [float]: latitude of satellite (deg)
        alts [float]: altitude of satellite (km)
        ephtimes [float]: Truncated JD (MET) times since 'era' (seconds)
        tx_lon [float]: Longitude of transmitter (degrees -180 to 180)
        tx_lat [float]: Latitude of transmitter (degrees -90 to 90)
        freq [float]: frequency of transmitted wave [Hz]
        omega_p [float]: plasma frequency of peak plasma in region [Rad/s]
        *REMOVED* omega_c [float]: cyclotron frequency of peak plasma in region [Rad/s]

    ***RETURNS***
        angles [float]: Aspect angles 
        nplus [float]: positive Appleton-Hartree solution (O-mode index of ref)
        nminus [float]: negative AH solution (X-mode index of ref)

    """
    kvecs,bvecs,angles = get_kb_angle(lons,lats,alts,ephtimes,tx_lon=tx_lon,tx_lat=tx_lat)
    #bvecs,kvecs,angles = get_kb_ottawa_angle(lons,lats,alts,ephtimes)
    angles_rad = np.deg2rad(angles)
    #TODO: allow user to input nu_e to this?
    X,Y,Z = appleton_coeffs(omega_c,omega_p,freq,0.0) 
    np_sq,nm_sq = appleton_hartree(X,Y,Z,angles_rad)
    nplus = np.sqrt(np_sq)
    nminus = np.sqrt(nm_sq)
    return angles_rad, nplus, nminus 

# ------------------------------------------------------------------------------
#                         Telemetry and Ephemeris Functions
#                         ---------------------------------
#       Functions for calculating speeds, directions, and angles during an ePOP
#       pass.
# ------------------------------------------------------------------------------

def get_kvecs(glon, glat, altitude, tx_lon=-75.552, tx_lat=45.403, tx_alt=0.07):
    """
    This function takes a satellite ephemeris point(s) as input, calculating
    the straight-line k-vector from a transmitter in terms of North (x), 
    East (y), Down (z) components.

    Default values for the transmitter parameters are the coordinates of the
    Ottawa transmitter from the 18-22 April 2016 RRI experiments.

    ***PARAMS***
        glon (float or float array): longitude(s) of final point(s) from transmitter
        glat (float or float array): latitude(s) of final point(s) from transmitter
        altitude (float or float array): altitude(s) of final point(s) from transmitter
        [tx_lon] (float): longitude of transmitter [deg]
        [tx_lat] (float): latitude of transmitter [deg]
        [tx_alt] (float): altitude of transmitter [km]

    ***RETURNS***
        kv (float or float array): the vector(s) from transmitter to input point(s) 

    WORKS WITH SINGLE POINTS AND NDARRAYS/LISTS
    """
    are_vectors = True if (type(glon)==list or type(glon)==np.ndarray) else False

    init_bearing,final_bearing = get_bearing(tx_lon, tx_lat, glon, glat) 

    # In spherical coordinates, subtracting vectors doesn't get us the path
    # from one point to another along the surface. To accomplish that, we use
    # bearings and elevation angles.
    
    init_bearing,final_bearing = get_bearing(tx_lon, tx_lat, glon, glat)
    if are_vectors:
        txlons = np.ones(len(glon))*tx_lon
        txlats = np.ones(len(glat))*tx_lat
        txalts = np.ones(len(altitude))*tx_alt
        elev_angle = get_elevation_angle(txlons, txlats, txalts, glon, glat, altitude)
    else:
        elev_angle = get_elevation_angle(tx_lon, tx_lat, tx_alt, glon, glat, altitude)

    # for calculation of vector components, we need theta (the bearing in the
    # N-E plane) at CASSIOPE, and phi, the angle-from-down (toward the N-E plane)
    theta = final_bearing
    phi = 90. + elev_angle
    kx = np.sin(np.deg2rad(phi))*np.cos(np.deg2rad(theta))
    ky = np.sin(np.deg2rad(phi))*np.sin(np.deg2rad(theta))
    kz = np.cos(np.deg2rad(phi)) 

    if are_vectors:
        kv = np.array([ (kx[i],ky[i],kz[i]) for i in range(len(kx))])
    else:
        kv = np.array((kx,ky,kz))
    return kv

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


    **With production of 'get_igrf' function which can handle array inputs, 
        this function is now R E D U N D A N T**

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


def get_igrf(lons,lats,alts,ephtimes):
    """
    This function uses the IGRF model in DavitPy to calculate the magnetic
    field vector at given near-Earth locations at a particular times.

    Currently, the function just takes geographic longitude and latitude. It 
    shouldn't be difficult to extend its functionality to accept magnetic
    coordinates as well though.

    **PARAMS**
    lons (Float ndarray): Geographic longitude of point
    lats (Float ndarray): Geographic latitude of point
    alts (Float ndarray): Altitude from Earth's surface in km
    ephtimes (Float ndarray): Time of interest (truncated julian time MET)
        (Time required because magnetic North's location varies year-to-year)

    """
    from davitpy.models import igrf
    from davitpy import utils
    itype = 1 #Geodetic coordinates
    stp = 1. 
    ifl = 0
    times = ephems_to_datetime(ephtimes)
    B_igrf = np.zeros((len(times),3))
    for i, time in enumerate(times):
        date = utils.dateToDecYear(time)
        lon = lons[i]
        lat = lats[i]
        alt = alts[i]
        xlti, xltf, xltd = lat, lat, stp
        xlni, xlnf, xlnd = lon, lon, stp 
        # Call fortran subroutine
        lat,lon,d,s,h,bx,by,bz,f = igrf.igrf11(itype,date,alt,ifl,xlti,xltf,xltd,xlni,xlnf,xlnd)
        B_igrf[i,:] = np.array((bx[0],by[0],bz[0]))
    return np.array(B_igrf)   

def get_aspect_angle(kv,bv,vector_inputs=False):
    """
    Given b-vector(s) and k-vector(s), computes angle between them (aspect angle).

    ***PARAMS***
        kv: single vector or list/array of vectors showing the propagation direction of wave
        bv: single vector or list/array of vectors for magnetic field
        [vector_inputs]: boolean flag indicating to treat vector inputs

    ***RETURNS***
        ang_deg: aspect angle [deg]
    """
    if vector_inputs:
        ang_deg = []
        for i,kv_i in enumerate(kv):
            bv_i = bv[i]
            prod_i = np.dot(kv_i,bv_i)/(np.linalg.norm(kv_i)*np.linalg.norm(bv_i))
            ang_deg_i = np.degrees(np.arccos(prod_i))
            ang_deg.append(ang_deg_i)
        ang_deg = np.array(ang_deg)
    else:
        prod = np.dot(kv,bv)/(np.linalg.norm(kv)*np.linalg.norm(bv)) 
        ang_deg = np.degrees(np.arccos(prod))
    return ang_deg

def get_plasma_intersection(lon,lat,alt,plasma_alt=300.,tx_lon=-75.552, tx_lat=45.403, tx_alt=0.07):
    """
    This function finds where a ray from a transmitter toward a satellite intersects the 
    peak plasma in the middle.

    ***PARAMS***
    Satellite ephemeris point(s): lon,lat,alt (deg, deg, km)
    Transmitter location [optionally]: tx_lon, tx_lat, tx_alt (deg, deg, km)
    Altitude of peak plasma density: plasma_alt (km.) 

    ***RETURNS***
    plasma_lon (float): longitude of plasma intersection(s)
    plasma_lat (float): latitude of plasma intersection(s) 

    """
    are_vectors = True if (type(lon)==list or type(lon)==np.ndarray) else False
    #lon = (lon + 360.) % 360.
    #tx_lon = (tx_lon + 360.) % 360.
    dist = haversine(lon,lat,tx_lon,tx_lat)
    if dist > 2500.:
        print("This approximation isn't valid for large distances")
        print("dist: ",dist)
        return (-1,-1)
    if plasma_alt > np.min(alt):
        print("Input altitudes are too low for the plasma")
        print('plasma_alt: ',plasma_alt)
        print('alt: ',alt) 
        return (-1,-1)
    if are_vectors:
        tx_lon = tx_lon*np.ones(len(lon))
        tx_lat = tx_lat*np.ones(len(lat))
        tx_alt = tx_alt*np.ones(len(alt))
    x = (plasma_alt/alt)*dist
    #print ('x, dist: ',x,dist)
    bearing,__ = get_bearing(tx_lon, tx_lat, lon, lat) #only need initial bearing
    #print ('bearing: ',bearing)
    delta_EW = x*np.sin(np.deg2rad(bearing)) 
    delta_lon = delta_EW*360./(2*np.pi*6371.*np.sin(np.deg2rad(lat))) # convert to longitude (deg)
    #delta_lon = delta_EW/(6371.*np.sin(np.deg2rad(lat))) # convert to longitude (deg)

    delta_NS = x*np.cos(np.deg2rad(bearing)) # small bearings=mostly northward-> mostly delta_lat
    delta_lat = delta_NS*360./(2*np.pi*6371.)
    #print('delta_EW,delta_NS: ',delta_EW,delta_NS)
    #print('delta_lon,delta_lat : ',delta_lon, delta_lat)

    plasma_lon = tx_lon + delta_lon
    plasma_lat = tx_lat + delta_lat
    #print('plasma_lon,plasma_lat: ',plasma_lon,plasma_lat)
    return (plasma_lon, plasma_lat)

def faraday_pass(lons,lats,alts,ephtimes,densities_arr,densities_lats,tx_lon=-75.552,tx_lat=45.503,tx_alt=0.07):
    """
    Performs faraday_trace() for each ephemeris point in a pass.
    """
    phase_integrals = []
    ql_integrals = []
    tecs = []
    mean_bcs = []
    for i,alt in enumerate(alts):
        lon = lons[i]
        lat = lats[i]
        ephtime = ephtimes[i]
        TEC,bcs,ql_integral,phase_integral = faraday_trace(lon,lat,alt,ephtime,densities_arr,densities_lats,tx_lon,tx_lat,tx_alt)
        tecs.append(TEC)
        mean_bcs.append(np.mean(bcs))
        ql_integrals.append(ql_integral)
        phase_integrals.append(phase_integral)
    return tecs,mean_bcs,ql_integrals,phase_integrals

def faraday_trace(lon,lat,alt,ephtime,densities_arr,densities_lats,tx_lon=-75.552,tx_lat=45.503,tx_alt=0.07,freq=1.0422E7):
    """
    Does a ray trace from a transmitter to a point to calculate faraday rotation.

    This ray trace also calculates individual terms along the way (line 
    integrals of electron density, B cos(theta), and their product).

    Directly computes indices of refraction for each voxel to allow a mostly 
    unapproximated calculation of phase.

    """
    print("Tracing from transmitter at (75W,45N,70m Alt) to ({0},{1},{2}km elevation)".format(lon,lat,alt)) 

    ang_deg = get_elevation_angle(tx_lon,tx_lat,tx_alt,lon,lat,alt)
    # path length through plasma voxel is r = y/sin(theta), y = 1E3 (1km)
    dist = 1E3/np.sin(np.deg2rad(ang_deg)) 
    kv = get_kvecs(lon,lat,alt)
    time = ephem_to_datetime(ephtime)

    # Plasma density profiles start at 60 km and go up to 559 km
    ndiffs = []
    bcs = []
    TEC = 0.
    phase_integral = 0.
    ql_integral = 0.

    pl_alts = np.arange(60.,alt)

    for pl_alt in pl_alts: # Go from 60 km altitude up to satellite altitude
        plon,plat = get_plasma_intersection(lon,lat,alt,plasma_alt=pl_alt,tx_lon=tx_lon,tx_lat=tx_lat,tx_alt=tx_alt)

        n_e = data_utils.get_density(plon, plat, pl_alt, densities_arr, densities_lats)
        TEC += dist*n_e

        omega_p = plasma_freq(n_e)
        bv = get_bvec(plon,plat,pl_alt,time)
        omega_c = cyclotron_freq(1E-9*np.linalg.norm(bv))

        X,Y,Z = appleton_coeffs(omega_p,omega_c,freq,0.)
        ang_deg = get_aspect_angle(kv,bv)

        b_cos_theta = np.abs(np.linalg.norm(bv)*np.cos(np.deg2rad(ang_deg)))
        bcs.append(b_cos_theta)
        ql_integral += dist*b_cos_theta*n_e

        nplus,nminus = appleton_hartree(X,Y,Z,np.deg2rad(ang_deg))
        phase_integral += dist*(nplus-nminus)

    factors = (2*np.pi*freq)/(2*3.00E8)
    phase_integral *= factors # omega/2c factor for the integral
    return TEC,bcs,ql_integral,phase_integral

def get_kb_angle(lons, lats, alts, ephtimes, tx_lon=-75.552, tx_lat=45.503, tx_alt=0.07):
    """
    This function takes the ephemeris data as arguments and using the Ottawa
    NRCAN transmitter, compares the k (line of sight) vector and the IGRF B 
    field vector to determine the relative angle between the propagating radio
    wave and the B-field in the ionospheric plasma.

    The purpose of this is to get an idea of which mode the radio wave is 
    propagating under (and if it's undergoing Faraday rotation, etc.)

    *** PARAMS ***
    lons (np.array of floats): ground-track geographic longitude
    lats (np.array of floats): ground-track geographic latitude
    alts (np.array of floats): altitude of CASSIOPE 
    ephtimes (np.array of floats): in truncated JD (MET), seconds since era
    [tx_lon] (np.array of floats): optional transmitter longitude
    [tx_lat] (np.array of floats): optional transmitter latitude
    [tx_alt] (np.array of floats): optional transmitter altitude

    *** RETURNS ***
    bvecs (np.array of float triplets): the IGRF B field at each ephemeris point
    kvecs (np.array of float triplets): the tx-to-rx(sat) vector at each ephemeris point
    angles (np.array of floats): the aspect angle (angle between respective bvecs and kvecs)
    """
    # Extension to allow single aspect angles to be calculated
    import numbers
    if isinstance(lons, numbers.Number):
        #print("Inputting single query")
        kv = get_kvecs(lons,lats,alts)
        bv = get_bvec(lons,lats,alts,ephem_to_datetime(ephtimes))
        prod = np.dot(kv,bv)/(np.linalg.norm(bv)*np.linalg.norm(kv))
        ang = np.rad2deg(np.arccos(prod))
        return bv,kv,ang 

    # Otherwise, for vector inputs:
    times = ephems_to_datetime(ephtimes)
    angles = []
    bvecs = []
    kvecs = []
    for i,lon in enumerate(lons):
        lat = lats[i]
        alt = alts[i]
        time = times[i]
        bvec = get_bvec(lon,lat,alt,time)
        kvec = get_kvecs(lon,lat,alt,tx_lon=tx_lon,tx_lat=tx_lat,tx_alt=tx_alt)
        #print "B vector: " + str(bvec)
        #print "K vector: " + str(kvec)
        # Take the dot product, divide out the magnitude of the vectors
        prod = np.dot(kvec,bvec)/(np.linalg.norm(bvec)*np.linalg.norm(kvec))
        angle = np.rad2deg(np.arccos(prod)) 
        angles.append(angle)
        bvecs.append(bvec)
        kvecs.append(kvec) 
    return np.array(bvecs),np.array(kvecs),np.array(angles)

def get_ramdirs(glon, glat, altitude):
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

def get_closest_approach(lons, lats, alts, tx_lon=OTTAWA_TX_LON, tx_lat=OTTAWA_TX_LAT, tx_alt=OTTAWA_TX_ELEV):
    """
    Uses a haversine formula-based distance calculation to find the closest 
    approach of a set of points describing a ground-track.

    *** PARAMS *** 
    glons (np.array of floats): longitudes (deg) of ground-track points
    glats (np.array of floats): latitudes (deg) of ground-track points
    alts (np.array of floats): altitudes (in km)

    ** RETURNS **
    index_shortest (integer): index of the lat/lon pair where the ottawa transmitter is closest
    dists (float): the calculated distances
    """

    # *** FINDING THE CLOSEST APPROACH ***
    # Using the Haversine formula (in a function in script_utils.py), the closest
    # approach is determined by brute force.
    dists = []
    longdists = []
    latdists = []
    dist_shortest = sys.maxint #Initially set this very high
    for i in range(np.size(lons)):
        dist_before_alt = haversine(lons[i], lats[i], tx_lon, tx_lat)
        delt_alt = alts[i] - tx_alt
        dist = np.linalg.norm((dist_before_alt,delt_alt))
        # Initially I took a quick and dirty approach to find the point of smallest
        # Euclidean distance in terms of latitudes and longitudes.
        #longdist = abs(geog_longs[i] - OTTAWA_TX_LON) # difference of longitudes
        #latdist = abs(geog_lats[i] - OTTAWA_TX_LAT) # difference of latitudes
        #dist = np.sqrt(longdist*longdist + latdist*latdist)  
        dists.append(dist)
        if dist < dist_shortest:
            dist_shortest = dist
            index_shortest = i
    return index_shortest, dists 

def get_kdip_angles(lons,lats,alts,ephtimes,pitch,yaw,roll):
    """
    Call this function to retrieve a list of the direction of the RRI dipole plane 
    in N-E-down coordinates for the first n-1 ephemeris points provided, as well
    as the relative angle between the K_los vector from the Ottawa transmitter
    and the direction of the RRI dipole plane.

    """
    #TODO: VALIDATE THIS. SEEMS LIKE WE'RE OFF BY 90 DEGREES

    # vs is a vector of ram directions in N, E, Down coordinates
    vs,dists = get_ramdirs(lons,lats,alts)
   
    # kvecs will be k_LOS vectors in N, E, Down coordinates 
    kvecs = []
    for i in range(lons.__len__()):
        kvec = get_kvecs(lons[i],lats[i],alts[i])
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
    # doesn't account for curvature of earth! Future goal
    if type(arcdist)==list or type(arcdist)==np.ndarray:
        if (arcdist > 2500.).any():
            print("**This approximation won't work for points this far apart**")
            return -1.
    elif arcdist > 2500.:
        print("**This approximation won't work for points this far apart**")
        return -1.

    elev_angle = np.rad2deg(np.arctan2(delta_alt,arcdist)) 
    return elev_angle

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

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    from math import radians
    # convert decimal degrees to radians 
    if type(lon1)==list or type(lon1)==np.ndarray:
        lon1 = np.deg2rad(lon1)
        lat1 = np.deg2rad(lat1)
        lon2 = np.deg2rad(lon2)
        lat2 = np.deg2rad(lat2)
    else:
        lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]]) 


if __name__=="__main__":

    """ TESTING INDEX OF REFRACTION-RELATED FUNCTIONS """
    freq10 = 1.0422E7   #[Hz]             
    freq12 = 1.2500E7   #[Hz]
    txlon = -75.552    #[Deg]
    txlat = 45.403     #[Deg]

    date_string="20160421" 
    filename, __ = data_utils.get_ottawa_data(date_string)
    omega_p = 2*np.pi*5.975E6 # [rad/s]
    fof2_alt = 260. # [km] 
    
    lons,lats,alts,ephtimes = data_utils.get_rri_ephemeris(filename)
    b_igrf = magnet_data.get_igrf(lons,lats,alts,ephtimes)

    # Test regular and improved cyclotron frequency functions
    '''
    I'm primarily checking that the functions don't break in general, but the
    cases below check particular numbers in order to detect whether the data or
    the processing unexpectedly changes after some updates to things.
    ''' 
    omega_c1 = improved_cyclotron_freq(lons,lats,alts,ephtimes,fof2_alt=280.)
    omega_c2 = cyclotron_freq(b_igrf)
    if np.round(np.log10(omega_c2)) != 7. or np.round(np.log10(omega_c1)) != 7.:
        print("Error with cyclotron frequency calculations")

    # Test plasma_freq(n_e) #TODO:

    # Test basic_ionosphere_params
    basic_omega_c, basic_omega_p, basic_l_d = basic_ionosphere_params(altitude=300.)
    if np.round(np.log10(basic_omega_p),1)!=7.8 or np.round(np.log10(basic_omega_c),1) != 7.0:
        print("Error with basic ionosphere params")

    # Test appleton_coeffs
    # X should be about 0.25-0.3, Y about 0.1, Z about 0
    X,Y,Z = appleton_coeffs(omega_c1, omega_p, freq10, 0.) 
    if np.round(X,1) != 0.33 or np.round(Y,2) != 0.13 or np.round(Z,2) != 0.:
        print("Error with appleton_coeffs()")

    # Test appleton_hartree
    theta = np.array( range(315) )/100.
    nplus1,nminus1 = appleton_hartree(X,Y,Z,theta)
    if np.round(abs(nplus1[0]),4) != 0.8416 or np.round(abs(nminus1[0]),4) != 0.7896:
        print("Error with appleton_hartree()")

    # Test txpass_to_indices 
    angles_r,nplus2,nminus2 = txpass_to_indices(lons,lats,alts,ephtimes,txlon,txlat,freq10,omega_p,omega_c1)
    if np.round(np.mean(angles_r),4) != 2.3198 or np.round(abs(nplus2[0]),4) != 0.9138:
        print("Error with txpass_to_indices()")



    """ TESTING DIRECTIONAL/EPHEMERIS-RELATED FUNCTIONS"""

    # Test get_bearing()
    # Correctness:
    tst_bearing1,__ = get_bearing(-106.,52.,-96.,52.)
    tst_bearing2,__ = get_bearing(-106.,52.,-106.,42.)
    if np.round(tst_bearing1,1) != 86.1 or np.round(tst_bearing2,1) != 180.:  
        print("Error calculating bearings")
    # Parallelization:
    a = np.array([100.,100.,100.])
    b = np.array([45.,45.,45.])
    c = np.array([35.,55.,-45.])
    bearings_i,bearings_f = get_bearing(a,b,a,c)
    if (bearings_i != [180.,0.,180.]).any(): # if any entries incorrect...
        print("Error with parallelized bearings calculation")

    # Test haversine:
    hav = haversine(50,50,60,60)
    if int(hav)!=1278: #corroborated with other calculators
        print("Error with haversine()")
    # check parallelization
    havs = haversine(a,b,a,c)
   
    # Test elevation_angle: 
    elev_a = get_elevation_angle(-106.,52.,0.,-96.,52.,300.)
    if np.round(elev_a,2) != 23.68:
        print("Error with get_elevation_angle()")

    # Test rotation_matrix:
    # For now I'm pretty confident from previous testing, let this be a #TODO

    # Test get_kvecs:
    kvs = get_kvecs(lons,lats,alts) # From Ottawa by default. 
    # I've already checked the validity of the April 21st data. This test relies on it:
    if (np.round((np.mean(kvs[:,0]),np.mean(kvs[:,1]),np.mean(kvs[:,2])),2) != [0.09,0.05,-0.65]).any():
        print("Error with get_kvecs()")

    # Test get_igrf(): already called get_igrf earlier in this if __name__ block
    #b_igrf = get_igrf(lons,lats,alts,ephtimes)
    if (np.round(b_igrf[0]) != np.array([17858.,-3148.,39417.])).any():
        print("Error with get_igrf()")

    # Test get_aspect_angle():
    kv_tst = np.array((1,2,3))
    bv_tst = np.array((2,30,7))

    # Test get_plasma_intersection():
    plasma_lon,plasma_lat = get_plasma_intersection(-75.552,48.403,370.) #directly north
    if np.round(plasma_lon,1)!=284.4 or np.round(plasma_lat,1)!=47.8:
        print("Error with get_plasma_intersection")

    # Test get_kb_angle():
    # Already tested successfully earlier by doing tx_pass_indices
    #bvs,kvs,angles_r = get_kb_angle

    # Test get_kdip_angles():
    #TODO: test and *validate*!
    lons,lats,alts,ephtimes,mlons,mlats,mlts,pitch,yaw,roll = data_utils.get_rri_ephemeris_full(filename)
    dipole_dirs, kdip_angles = get_kdip_angles(lons,lats,alts,ephtimes,pitch,yaw,roll)
 
    # Test get_ramdirs():
    vs,dists = get_ramdirs(lons,lats,alts) 
    if (np.round(vs[0][0],3)!=7.442) or (np.round(dists[-1],3)!=7.419):
        print("Error with get_ramdirs")

    # Test get_closest_approach:
    index_closest,dists = get_closest_approach(lons,lats,alts)
    if (index_closest!=103) or np.round(dists[-1],2)!=1049.97:
        print("Error with get_closest_approach")

    # Test faraday_trace
    datarr,datlats = data_utils.load_density_profile('./data/20160418-densities.txt')
    tec,bcs,__,__ = faraday_trace(lons[0],lats[0],alts[0],ephtimes[0],datarr,datlats)
    if np.round(np.log10(tec),3)!=17.198: # validated calculations by noting similarity of numbers with Rob's
        print("Error with faraday_trace") 



    # Some extra stuff I want to do for analysis at the bottom here...
    Roblons,Roblats,Robalts = data_utils.load_rob_ephemeris('./data/satcoords_20160418.txt')
    lons18,lats18,alts18,ephtimes18 = data_utils.get_rri_ephemeris(data_utils.get_ottawa_data('20160418')[0])    

    plons=[]
    plats=[]
    i = 19 
    lon18 = lons18[i]
    lat18 = lats18[i]
    alt18 = alts18[i]
    ephtime18 = ephtimes18[i]
    plalts = np.arange(60.,alt18)
    for plalt in plalts:                                                   
        plon,plat = get_plasma_intersection(lon18,lat18,alt18,plasma_alt=plalt)
        plons.append(plon)
        plats.append(plat)

    #tecs,mean_bcs,ql_integrals,phase_integrals = faraday_pass(lons18,lats18,alts18,ephemtimes18,datarr,datlats)
    # plt.plot(phase_integrals); plt.show()

    print("Tests complete!")
