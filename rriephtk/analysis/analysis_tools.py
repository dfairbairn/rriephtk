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
import logging

import __init__ # 
import rriephtk.utils.data_utils as data_utils
from rriephtk.utils.data_utils import ephem_to_datetime, ephems_to_datetime 

import magnet_data

OTTAWA_TX_LON = -75.552
OTTAWA_TX_LAT = 45.403
OTTAWA_TX_ELEV = 0.070 # Ottawa 70m elevation - small compared to satellite

# Physical constants:
EL_CHARGE = 1.602E-19   #[C]
EL_MASS = 9.109E-31 #[kg]
EPS0 = 8.8542E-12   #[A*s/(V*m)]
EARTH_RAD = 6371.   #[km]

logging.basicConfig(filename='analysis-tools.log',level=logging.WARNING)

# -----------------------------------------------------------------------------
#                       Index of Refraction-related Code
#                       --------------------------------
#    The purpose of this code is the calculation of index of refraction vs
#    incidence angle for radio wave propagation direction compared to the
#    ionospheric plasma's orientation/local B field
# -----------------------------------------------------------------------------

def basic_ionosphere_params(altitude=300.):
    """ Returns cyclotron frequency, debye length, and plasma frequency in
    ionosphere for a given altitude.

    *PARAMS*:
        [altitude]: In km. Must be in ionosphere for this model to work.
                    Default: 300 km.

    ANGULAR FREQUENCIES BTW

    *RETURNS*:
        omega_p: Plasma frequency (angular)
        omega_c: Cyclotron frequency (angular)
        l_d: Debye length

    Based on:
    https://smallsats.org/2013/03/25/debye-length-plasma-frequency-and-cyclotron-frequency/
    """

    if altitude > 400. or altitude < 200.:
        logging.error("Bad altitude parameter")
        return None

    # Assumed parameters:
    me   = 9.109E-31            #[kg]           Electron rest mass
    mp   = 1.673E-27            #[kg]           Proton rest mass
    eps0 = 8.8542E-12           #[A*s/(V*m)]    Permittivity
    e    = 1.602E-19            #[C]            Elementary charge
    Re   = 6.37E6               #[m]            Earth radius
    Md   = 8E15                 #[T*m3]         Earth magnetic dipole moment

    n   = 1E12                     #[m-3]  Electron density < HIGH >
    kTe = 0.1                      #[eV]   Thermal energy
    r   = Re + altitude*1000.      #[km]   Radial dist from Earth’s center
    omega_p  = ((e**2)*n/(eps0*me))**0.5  #[Hz]   Plasma frequency
    l_d  = (eps0*kTe/(n*e))**0.5     #[m]    Debye length
    B   = 2*Md/(r**3)                #[T]    Earth’s magnetic field
    omega_c  = e*B/(me)            #[Hz]   Cyclotron frequency

    return omega_c, omega_p, l_d

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

    Basically either uses a single B field for a cyclotron frequency, or takes
    the average of an array and uses that to find a cyclotron frequency.

    *** PARAMS ***
        B [float] or [list/ndarray of floats]: magnetic field. If in list form,
            this function takes the average magnitude of B's [T]
         *Note that B is expected to be in tesla (e.g. around 50000nT is 50000.E-9)

    ***RETURNS***
        omega_c [float]: cyclotron frequency (rad/s)
    """
    e    = 1.602E-19            #[C]            Elementary charge
    me   = 9.109E-31            #[kg]           Electron rest mass
    if type(B)==np.ndarray or type(B)==list:
        # Convert with 1E-9 for [m^3] vs [km^3]
        B_mag = 1E-9*np.mean([np.linalg.norm(j) for j in B ])
    elif isinstance(B, float):
        B_mag = B
    else:
        logging.error("Unexpected data type for B field parameter")
        logging.error("{0}\t{1}".format(B,type(B)))
        return None
    omega_c  = e*B_mag/(me)            #[Hz]   Cyclotron frequency
    return omega_c

def improved_cyclotron_freq(lons, lats, alts, ephtimes, fof2_alt=250.):
    """
    Takes CASSIOPE ephemeris (lons,lats,alts,ephtimes) and uses IGRF to
    calculate predicted B-field strength, returning a mean cyclotron freq
    based on the mean B field over the course of the pass.

    *** PARAMS ***
    lons: longitude (deg)
    lats: latitude (deg)
    alts: altitude (km)
    ephtimes: times 
    [fof2_alt] (default: 250.): height of the fof2 peak in km. Subroutine uses 
                            the fof2 peak alt. as the location at which to get
                            model's B field
    
    *** RETURNS ***
    omega_c: average cyclotron frequency encountered by incident radio waves in
            ionosphere (in the sliver of ionosphere where lots of stuff happens)
    """
    import magnet_data
    e =  1.602E-19            #[C]            Elementary charge
    me = 9.109E-31            #[kg]           Electron rest mass

    # attempt to find B field within peak plasma density region at fof2 alt
    B = magnet_data.get_igrf(lons, lats, fof2_alt*np.ones(len(alts)), ephtimes)
    # Convert with 1E-9 for [m^3] vs [km^3]
    B_mag = 1E-9*np.mean([np.linalg.norm(j) for j in B ])
    omega_c = e*B_mag/me
    return omega_c


def appleton_coeffs(omega_c, omega_p, freq, nu_e):
    """
    Gets the Appleton Hartree X, Y, Z coefficients given the ionospheric params.

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
    return X, Y, Z


def appleton_hartree(X, Y, Z, theta):
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
    
    NOTE: CURRENTLY USES APPROXIMATIONS ASSUMING FLAT EARTH (ONLY
    REASONABLE NEAR THE TRANSMITTER). CAN BE OFF BY UP TO 8 DEGREES
    FOR THE APRIL 2016 PASSES.
    
    USE 'GET_KVEC2', 'GET_BVEC2' FOR IMPROVED PERFORMANCE.
    """

    term1 = 1. - 1j*Z
    term2 = ((Y*np.sin(theta))**2)/(2*(1 - X - 1j*Z))
    term3 = np.sqrt(((Y*np.sin(theta))**4)/(4*((1 - X - 1j*Z)**2)) +
            (Y*np.cos(theta))**2)
    nplus_squared = 1 - X/(term1 - term2 + term3)
    nminus_squared = 1 - X/(term1 - term2 - term3)
    nplus = np.sqrt(nplus_squared)
    nminus = np.sqrt(nminus_squared)
    return nplus, nminus

def txpass_to_indices(lons, lats, alts, ephtimes,
                      tx_lon, tx_lat, freq, omega_p, omega_c, nu_e=0.0):
    """
    Calculates the expected indices of refraction at the fof2 peak for radio
    waves directed from a specified transmitter toward CASSIOPE.

    Each satellite ephemeris point ([lon, lat, alt, ephtime]) receives aspect
    angles (from get_kb_angle) and corresponding indices of refraction for the
    ionospheric plasma at the given frequency.

    *** PARAMS ***
        lons [float]: longitude of satellite (deg)
        lats [float]: latitude of satellite (deg)
        alts [float]: altitude of satellite (km)
        ephtimes [float]: Truncated JD (MET) times since 'era' (seconds)
        tx_lon [float]: Longitude of transmitter (degrees -180 to 180)
        tx_lat [float]: Latitude of transmitter (degrees -90 to 90)
        freq [float]: frequency of transmitted wave [Hz]
        omega_p [float]: plasma frequency of peak plasma in region [Rad/s]
        omega_c [float]: cyclotron frequency of peak plasma in region [Rad/s]

    ***RETURNS***
        angles [float]: Aspect angles
        nplus [float]: positive Appleton-Hartree solution (O-mode index of ref)
        nminus [float]: negative AH solution (X-mode index of ref)

    """
    kvecs, bvecs, angles = get_kb_angle(lons, lats, alts, ephtimes,
                                      tx_lon=tx_lon, tx_lat=tx_lat)
    #bvecs, kvecs, angles = get_kb_ottawa_angle(lons, lats, alts, ephtimes)
    angles_rad = np.deg2rad(angles)
    X, Y, Z = appleton_coeffs(omega_c, omega_p, freq, nu_e)
    np_sq, nm_sq = appleton_hartree(X, Y, Z, angles_rad)
    nplus = np.sqrt(np_sq)
    nminus = np.sqrt(nm_sq)
    return angles_rad, nplus, nminus

# -----------------------------------------------------------------------------
#                         Telemetry and Ephemeris Functions
#                         ---------------------------------
#       Functions for calculating speeds, directions, and angles during an ePOP
#       pass.
# -----------------------------------------------------------------------------

def get_kvecs(lon, lat, altitude, tx_lon=-75.552, tx_lat=45.403, tx_alt=.07):
    """
    Uses North-East-Down coordinates to get a direction vector from a 
    transmitter to an ephemeris point in North-East-Down coordinates at
    the ephemeris point location.

    Default values for the transmitter parameters are the coordinates of the
    Ottawa transmitter from the 18-22 April 2016 RRI experiments.

    *** PARAMS ***
        lon (float or float array): longitude of point(s) from transmitter
        lat (float or float array): latitude of point(s) from transmitter
        altitude (float or float array): altitude of point(s) from transmitter
        [tx_lon] (float): longitude of transmitter [deg]
        [tx_lat] (float): latitude of transmitter [deg]
        [tx_alt] (float): altitude of transmitter [km]

    ***RETURNS***
        kv (float or float array): vector(s) from transmitter to input point(s)

    WORKS WITH SINGLE POINTS AND NDARRAYS/LISTS
    """
    vec_inp = True if (type(lon)==list or type(lon)==np.ndarray) else False

    init_bearing, final_bearing = get_bearing(tx_lon, tx_lat, lon, lat)

    # In spherical coordinates, subtracting vectors doesn't get us the path
    # from one point to another along the surface. To accomplish that, we use
    # bearings and elevation angles.

    init_bearing, final_bearing = get_bearing(tx_lon, tx_lat, lon, lat)
    if vec_inp:
        txlons = np.ones(len(lon))*tx_lon
        txlats = np.ones(len(lat))*tx_lat
        txalts = np.ones(len(altitude))*tx_alt
        elev_angle = get_elevation_angle(txlons, txlats, txalts,
                                         lon, lat, altitude)
    else:
        elev_angle = get_elevation_angle(tx_lon, tx_lat, tx_alt,
                                         lon, lat, altitude)

    # for calculation of vector components, we need theta (the bearing in the
    # N-E plane) at CASSIOPE, & phi, the angle-from-down (toward the N-E plane)
    theta = final_bearing
    phi = 90. + elev_angle
    kx = np.sin(np.deg2rad(phi))*np.cos(np.deg2rad(theta))
    ky = np.sin(np.deg2rad(phi))*np.sin(np.deg2rad(theta))
    kz = np.cos(np.deg2rad(phi))

    if vec_inp:
        kv = np.array([(kx[i], ky[i], kz[i]) for i in range(len(kx))])
    else:
        kv = np.array((kx, ky, kz))
    return kv

def get_kvec2(lon, lat, alt, ephtimes=None, 
              tx_lon=-75.552, tx_lat=45.403, tx_alt=.07):
    """
    Uses cartesian coordinates (GEO XYZ) to get absolute direction 
    vector between two (lon, lat, alt) points.

    *** PARAMS ***
        lon (float or array): longitude of point(s) from transmitter
        lat (float or array): latitude of point(s) from transmitter
        altitude (float or array): altitude of point(s) from transmitter
        [tx_lon] (float): longitude of transmitter [deg]
        [tx_lat] (float): latitude of transmitter [deg]
        [tx_alt] (float): altitude of transmitter [km]
    
    *** RETURNS ***
        kv (float or array): vector(s) from transmitter to input point(s)
       
    """
    import spacepy.coordinates as coord
    import spacepy.time as tm
    import datetime as dt
    b = coord.Coords([[alt + EARTH_RAD, lat, lon]],'GEO','sph')
    a = coord.Coords([[tx_alt + EARTH_RAD, tx_lat, tx_lon]], 'GEO', 'sph')
    if ephtimes is None:
        b.ticks = tm.Ticktock(dt.datetime.now()) # because it doesnt matter
        a.ticks = tm.Ticktock(dt.datetime.now()) # because it doesnt matter
    else:
        times = ephems_to_datetime(ephtimes)
        b.ticks = tm.Ticktock(times)
        a.ticks = tm.Ticktock(times)

    b = b.convert('GEO','car')
    a = a.convert('GEO','car')
    kv = (b.data - a.data)[0] 
    kv = kv/np.linalg.norm(kv)
    logging.info("get_kvecs2 result: ", kv)
    return kv

def get_bvec(lon, lat, altitude, time):
    """
    This function uses the IGRF model in DavitPy to calculate the B
    field vector at a given near-Earth location at a particular time.

    The field vector returned is in North-East-Down coordinates from
    the particular location's reference frame.

    ** PARAMS **
        lon (Float): Geographic longitude of point
        lat (Float): Geographic latitude of point
        alt (Float): Altitude from Earth's surface in km
        time (Datetime): Time of interest (magnetic North's location
                         varies year-to-year)

    *** RETURNS ***
        bv_ned (np.ndarray triplet of floats): magnetic field vector in
                the geographic North-East-Down coordinate system.

    **With production of 'get_igrf' function which handles array inputs,
        this is now R E D U N D A N T except for individual use**

    """
    from davitpy.models import igrf
    from davitpy import utils
    itype = 1 #Geodetic coordinates
    date = utils.dateToDecYear(time) # decimal year
    alt = altitude
    stp = 1. #
    # The IGRF function takes a grid of latitude and longitude points for which
    # to calculate the field. We just want one point.
    xlti, xltf, xltd = lat, lat, stp
    xlni, xlnf, xlnd = lon, lon, stp
    ifl = 0 # Main field
    # Call fortran subroutine
    lat_o, lon_o, d, s, h, bx, by, bz, f = igrf.igrf11(itype, date, alt, ifl,
                                           xlti, xltf, xltd, xlni, xlnf, xlnd)
    bv_ned = np.array((bx[0], by[0], bz[0]))
    return bv_ned

def get_bvec2(lon, lat, alt, time):
    """
    A method for getting the Earth's B-vector in GEO X-Y-Z coordinates.

    Uses the regular IGRF function for acquiring magnetic field, then 
    converts it using dir_ned2geo().

    *** PARAMS ***
        lon (Float): Geographic longitude of point
        lat (Float): Geographic latitude of point
        alt (Float): Altitude from Earth's surface in km
        time (Datetime): Time of interest (magnetic North's location
                         varies year-to-year)

    *** RETURNS ***
        bv_geo (np.ndarray triplet of floats): magnetic field vector in
                the geographic X Y Z coordinate system.
    """
    bvec = get_bvec(lon, lat, alt, time)
    bv_geo = dir_ned2geo((alt, lat, lon), bvec, time=time)
    logging.info('alt, lat, lon: ', alt, lat, lon)
    logging.info('B_NED = {ned},\tB_XYZ = {xyz}'.format(ned=bvec, xyz=bv_geo))
    return bv_geo

def get_igrf(lons, lats, alts, ephtimes):
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
        lat, lon, d, s, h, bx, by, bz, f = igrf.igrf11(
                itype, date, alt, ifl, xlti, xltf, xltd, xlni, xlnf, xlnd)
        # Save vector components to the B_igrf array
        B_igrf[i,:] = np.array((bx[0], by[0], bz[0]))
    return np.array(B_igrf)

def get_aspect_angle(kv, bv, vec_inp=False):
    """
    Given b-vec(s) and k-vec(s), computes angle between them (aspect angle).

    *** PARAMS ***
        kv: single vector or list/array of vectors showing the
            propagation direction of wave
        bv: single vector or list/array of vectors for magnetic field
        [vec_inp]: boolean flag indicating to treat vector inputs

    ***RETURNS***
        ang_deg: aspect angle [deg]
    """
    if vec_inp:
        ang_deg = []
        for i, kv_i in enumerate(kv):
            bv_i = bv[i]
            p_i = np.dot(kv_i, bv_i)/(np.linalg.norm(kv_i)*np.linalg.norm(bv_i))
            ang_deg_i = np.degrees(np.arccos(p_i))
            ang_deg.append(ang_deg_i)
        ang_deg = np.array(ang_deg)
    else:
        p = np.dot(kv, bv)/(np.linalg.norm(kv)*np.linalg.norm(bv))
        ang_deg = np.degrees(np.arccos(p))
    return ang_deg

def get_plasma_intersection(lon, lat, alt, plasma_alt=300., tx_lon=-75.552,
                            tx_lat=45.403, tx_alt=0.07):
    """
    This function finds where a ray from a transmitter toward a satellite
    intersects the peak plasma in the middle.

    *** PARAMS ***
    Satellite ephemeris point(s): lon, lat, alt (deg, deg, km)
    Transmitter location [optionally]: tx_lon, tx_lat, tx_alt (deg, deg, km)
    Altitude of peak plasma density: plasma_alt (km.)

    ***RETURNS***
    plasma_lon (float): longitude of plasma intersection(s)
    plasma_lat (float): latitude of plasma intersection(s)

    """
    vec_inp = True if (type(lon)==list or type(lon)==np.ndarray) else False
    #lon = (lon + 360.) % 360.
    #tx_lon = (tx_lon + 360.) % 360.
    dist = haversine(lon, lat, tx_lon, tx_lat)
    if dist > 2500.:
        logging.error("This approximation isn't valid for large distances")
        logging.error("dist: {0}".format(dist))
        return (-1,-1)
    if plasma_alt > np.min(alt):
        logging.error("Input altitudes are too low for the plasma")
        logging.error('plasma_alt: {0}'.format(plasma_alt))
        logging.error('alt: {0}'.format(alt))
        return (-1,-1)
    if vec_inp:
        tx_lon = tx_lon*np.ones(len(lon))
        tx_lat = tx_lat*np.ones(len(lat))
        tx_alt = tx_alt*np.ones(len(alt))
    x = (plasma_alt/alt)*dist
    #only need initial bearing
    bearing,__ = get_bearing(tx_lon, tx_lat, lon, lat)
    delta_EW = x*np.sin(np.deg2rad(bearing))
    delta_NS = x*np.cos(np.deg2rad(bearing))
    # convert to longitude (deg):
    delta_lon = delta_EW*360./(2*np.pi*6371.*np.sin(np.deg2rad(lat)))
    delta_lat = delta_NS*360./(2*np.pi*6371.)

    plasma_lon = tx_lon + delta_lon
    plasma_lat = tx_lat + delta_lat
    logging.info('delta_EW, delta_NS: {0},{1}'.format(delta_EW, delta_NS))
    logging.info('delta_lon, delta_lat: {0},{1}'.format(delta_lon, delta_lat))
    logging.info('plasma_lon, plasma_lat: {0},{1}'.format(plasma_lon, plasma_lat))
    return (plasma_lon, plasma_lat)

def faraday_pass(lons, lats, alts, ephtimes, densities_arr, densities_lats,
                 tx_lon=-75.552, tx_lat=45.503, tx_alt=0.07):
    """
    Performs faraday_trace() for each ephemeris point in a pass.

    *** PARAMS ***



    ***RETURNS***
        tecs
        mean_bcs [np.ndarray[float]]:
        ql_integrals [np.ndarray[float]]: integrals of B cos(theta)*N_e
        faraday_integrals [np.ndarray[float]]: integrals of difference
                                               in indices of refraction.
    """
    faraday_integrals = []
    ql_integrals = []
    tecs = []
    mean_bcs = []
    for i, alt in enumerate(alts):
        lon = lons[i]
        lat = lats[i]
        ephtime = ephtimes[i]
        TEC, bcs, ql_integral, faraday_integral = faraday_trace(
            lon, lat, alt, ephtime, densities_arr, densities_lats,
            tx_lon, tx_lat, tx_alt)
        tecs.append(TEC)
        mean_bcs.append(np.mean(bcs))
        ql_integrals.append(ql_integral)
        faraday_integrals.append(faraday_integral)
    return tecs, mean_bcs, ql_integrals, faraday_integrals

def faraday_trace(lon, lat, alt, ephtime, densities_arr, densities_lats,
                  tx_lon=-75.552, tx_lat=45.503, tx_alt=0.07, freq=1.0422E7):
    """
    Does a ray trace from a transmitter to a input point to calculate
    faraday rotation.

    This ray trace also calculates individual terms along the way (line
    integrals of electron density, B cos(theta), and their product).

    Directly computes indices of refraction for each voxel to allow a
    mostly unapproximated calculation of phase.

    """
    logging.info("Tracing from transmitter at (75W,45N,70m Alt) "
          "to ({0},{1},{2}km elevation)".format(lon, lat, alt))

    ang_deg = get_elevation_angle(tx_lon, tx_lat, tx_alt, lon, lat, alt)
    # path length through plasma voxel is r = y/sin(theta), y = 1E3 (1km)
    dist = 1E3/np.sin(np.deg2rad(ang_deg))
    kv = get_kvecs(lon, lat, alt)
    time = ephem_to_datetime(ephtime)

    # Plasma density profiles start at 60 km and go up to 559 km
    ndiffs = []
    bcs = []
    TEC = 0.
    phase_integral = 0.
    ql_integral = 0.

    pl_alts = np.arange(60., alt)

    for pl_alt in pl_alts: # Go from 60 km altitude up to satellite altitude
        plon, plat = get_plasma_intersection(lon, lat, alt, plasma_alt=pl_alt,
                                            tx_lon=tx_lon, tx_lat=tx_lat,
                                            tx_alt=tx_alt)

        n_e = data_utils.get_density(plon, plat, pl_alt, densities_arr,
                                     densities_lats)
        TEC += dist*n_e
        omega_p = plasma_freq(n_e)
        bv = get_bvec(plon, plat, pl_alt, time)
        # 1E-9 factor necessary for nT -> T
        omega_c = cyclotron_freq(1E-9*np.linalg.norm(bv))

        X, Y, Z = appleton_coeffs(omega_p, omega_c, freq,0.)
        ang_deg = get_aspect_angle(kv, bv)
        b_cos_theta = np.abs(np.linalg.norm(bv)*np.cos(np.deg2rad(ang_deg)))
        bcs.append(b_cos_theta)
        # 1E-9 factor necessary for nT -> T
        ql_integral += dist*b_cos_theta*n_e*1E-9

        nplus, nminus = appleton_hartree(X, Y, Z, np.deg2rad(ang_deg))
        phase_integral += dist*(nplus-nminus)

    factors1 = ((1.602E-19)**3)/(2*3.00E8*8.85E-12*(9.11E-31*2*np.pi*freq)**2)
    factors2 = (2*np.pi*freq)/(2*3.00E8)
    ql_integral *= factors1
    phase_integral *= factors2 # omega/2c factor for the integral
    return TEC, bcs, ql_integral, phase_integral

def get_kb_angle(lons, lats, alts, ephtimes,
                 tx_lon=-75.552, tx_lat=45.503, tx_alt=0.07):
    """
    This function takes the ephemeris data as arguments and using the
    Ottawa NRCAN transmitter, compares the k (line of sight) vector
    and the IGRF B field vector to determine the relative angle between
    the propagating radio wave and the Bfield in the ionospheric plasma.

    The purpose of this is to get an idea of which mode the radio wave
    is propagating under (and if it's undergoing Faraday rotation, etc.)


    *** PARAMS ***
      lons (np.array[float]): ground-track geographic longitude
      lats (np.array[float]): ground-track geographic latitude
      alts (np.array[float]): altitude of CASSIOPE
      ephtimes (np.array[float]): in truncated JD (MET),
                                  seconds since era
      [tx_lon] (np.array[float]): optional transmitter longitude
      [tx_lat] (np.array[float]): optional transmitter latitude
      [tx_alt] (np.array[float]): optional transmitter altitude

    *** RETURNS ***
      bvecs [np.array[(float, float, float)]: the IGRF B field at each
                                            ephemeris point
      kvecs [np.array[(float, float, float)]: the tx-to-rx(sat) vector
                                            at each ephemeris point
      angles [np.array[float]: the aspect angle (angle between
                               respective bvecs and kvecs)

    NOTE: CURRENTLY USES APPROXIMATIONS ASSUMING FLAT EARTH (ONLY
    REASONABLE NEAR THE TRANSMITTER). CAN BE OFF BY UP TO 8 DEGREES
    FOR THE APRIL 2016 PASSES.
    
    USE 'GET_KVEC2', 'GET_BVEC2' FOR IMPROVED ACCURACY.
    """

    # Extension to allow single aspect angles to be calculated
    import numbers
    if isinstance(lons, numbers.Number):
        logging.info("Inputting single query")
        kv = get_kvecs(lons, lats, alts)
        bv = get_bvec(lons, lats, alts, ephem_to_datetime(ephtimes))
        prod = np.dot(kv, bv)/(np.linalg.norm(bv)*np.linalg.norm(kv))
        ang = np.rad2deg(np.arccos(prod))
        return bv, kv, ang

    # Otherwise, for vector inputs:
    times = ephems_to_datetime(ephtimes)
    angles = []
    bvecs = []
    kvecs = []
    for i, lon in enumerate(lons):
        lat = lats[i]
        alt = alts[i]
        time = times[i]
        bvec = get_bvec(lon, lat, alt, time)
        kvec = get_kvecs(lon, lat, alt, 
                         tx_lon=tx_lon, tx_lat=tx_lat, tx_alt=tx_alt)
        logging.info("get_kb_angle B vector: {0}".format(bvec))
        logging.info("get_kb_angle K vector: {0}".format(kvec))
        # Take the dot product, divide out the magnitude of the vectors
        prod = np.dot(kvec, bvec)/(np.linalg.norm(bvec)*np.linalg.norm(kvec))
        angle = np.rad2deg(np.arccos(prod))
        angles.append(angle)
        bvecs.append(bvec)
        kvecs.append(kvec)
    return np.array(bvecs), np.array(kvecs), np.array(angles)

def get_ramdirs(glon, glat, altitude):
    """
    Computes the velocity components of the satellite based on its
    ephemeris points.

    *** PARAMS ***
    glon (np.array of floats): the longitudes of the satellite (deg)
    glat (np.array of floats): the latitudes of the satellite (deg)
    altitude (np.array of floats): the altitudes of the satellite (deg)
    time (doesnt even matter right now): The times (currently unused)

    *** RETURNS ***
    v (list of floats): velocity vector components (in NED) [km/s]
    dists (list of floats): full distances traveled in 1 sec in km

    *** VALIDATION ***
    Calculated velocity values agreed with velocity values from other
    sources

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
        # if lon1, lat1 -- lon2, lat2 is the hypotenuse, this is the vertical
        latdist = haversine(lon1, lat1, lon1, lat2)
        # and this is the horizontal
        londist = haversine(lon1, lat1, lon2, lat1)
        altdist = alt2 - alt1
        dist = np.linalg.norm((latdist, londist, altdist))
        v.append((latdist, londist, altdist))
        dists.append(dist)
    return v, dists

def get_closest_approach(lons, lats, alts, tx_lon=OTTAWA_TX_LON,
                         tx_lat=OTTAWA_TX_LAT, tx_alt=OTTAWA_TX_ELEV):
    """
    Uses a haversine formula-based distance calculation to find the
    closest approach of a set of points describing a ground-track.

    *** PARAMS ***
    glons (np.array of floats): longitudes (deg) of ground-track points
    glats (np.array of floats): latitudes (deg) of ground-track points
    alts (np.array of floats): altitudes (in km)

    ** RETURNS **
    index_shortest (integer): index of the lat/lon pair where the
                              ottawa transmitter is closest
    dists (float): the calculated distances
    """

    # *** FINDING THE CLOSEST APPROACH ***
    # Using the Haversine formula (in a function in script_utils.py), the
    # closest approach is determined by brute force.
    dists = []
    longdists = []
    latdists = []
    dist_shortest = sys.maxint #Initially set this very high
    for i in range(np.size(lons)):
        dist_before_alt = haversine(lons[i], lats[i], tx_lon, tx_lat)
        delt_alt = alts[i] - tx_alt
        dist = np.linalg.norm((dist_before_alt, delt_alt))
        # Initially, I took a quick and dirty approach to find the
        # point of smallest Euclidean distance in terms of latitudes
        # and longitudes.
        #longdist = abs(geog_longs[i] - OTTAWA_TX_LON) # difference of lons
        #latdist = abs(geog_lats[i] - OTTAWA_TX_LAT) # difference of lats
        #dist = np.sqrt(longdist*longdist + latdist*latdist)
        dists.append(dist)
        if dist < dist_shortest:
            dist_shortest = dist
            index_shortest = i
    return index_shortest, dists

def get_kdip_angles(lons, lats, alts, ephtimes, pitch, yaw, roll):
    """
    Call this function to retrieve a list of the directions of the RRI
    dipole plane in N-E-down coordinates for the first n-1 ephemeris
    points provided, as well as the relative angle between the K_los
    vector from the Ottawa transmitter and the direction of the RRI
    dipole plane.

    *** PARAMS ***
      lons (np.array[float]): ground-track geographic longitude
      lats (np.array[float]): ground-track geographic latitude
      alts (np.array[float]): altitude of CASSIOPE
      ephtimes (np.array[float]): in truncated JD (MET),
                                  seconds since era
      pitch (np.array[float]): 'up/down' spacecraft pitch from default
      yaw (np.array[float]): 'left/right' spacecraft yaw from default
      roll (np.array[float]): 'CW/CCW' rotational offset of spacecraft

    *** RETURNS ***
        dipole_dirs (np.array[(float, float, float)]: array of directions
            of the spacecraft at various points in N-E-D coords
        kdip_angles (np.array[float]): array of angles between the line
            of sight vector and the dipole_dirs

    """
    #TODO: VALIDATE THIS. SEEMS LIKE WE'RE OFF BY 90 DEGREES

    # vs is a vector of ram directions in N, E, Down coordinates
    vs, dists = get_ramdirs(lons, lats, alts)

    # kvecs will be k_LOS vectors in N, E, Down coordinates
    kvecs = []
    for i in range(lons.__len__()):
        kvec = get_kvecs(lons[i], lats[i], alts[i])
        kvecs.append(kvec)

    # Word of Gareth:
    # x is ram direction, z is nadir direction, y is Z cross X.
    xdirs = vs
    zdirs = np.array([(0,0,1) for i in range(xdirs.__len__())])
    ydirs = np.cross(zdirs, xdirs)

    # yaw: rot around z, pitch: rot around y, roll: rot around x
    # Assuming as seems to be confirmed in documentation that the
    # Dipole is in the x-direction on CASSIOPE, and thus, in default
    # position, towards the ram direction.
    dipole_dirs = []
    kdip_angles = []
    for i in range(xdirs.__len__()):
        yaw_rot = rotation_matrix(zdirs[i], np.deg2rad(yaw[i]))
        pitch_rot = rotation_matrix(ydirs[i], np.deg2rad(pitch[i]))
        roll_rot = rotation_matrix(xdirs[i], np.deg2rad(roll[i]))
        initial_dipole_vec = xdirs[i]
        intermed1 = np.dot(yaw_rot, initial_dipole_vec)
        intermed2 = np.dot(pitch_rot, intermed1)
        # After all the rotations, we have dip_dir
        dip_dir = np.dot(roll_rot, intermed2)
        dipole_dirs.append(dip_dir)
        # Determine what the relative angle is between k_LOS and dip_dir
        kdip_angle = np.arccos(np.dot(dip_dir, kvecs[i])/(
            np.linalg.norm(dip_dir)*np.linalg.norm(kvecs[i])))
        kdip_angles.append(np.rad2deg(kdip_angle))
    return dipole_dirs, kdip_angles

def get_elevation_angle(lon1, lat1, alt1, lon2, lat2, alt2):
    """
    Takes two lon/lat/alt points and determines the elevation angle
    necessary for getting from the first point to the second.

    *** PARAMS ***
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
            logging.error("*Approximation won't work for this distance*")
            return -1.
    elif arcdist > 2500.:
        logging.error("**Approximation won't work for points this distance**")
        return -1.

    elev_angle = np.rad2deg(np.arctan2(delta_alt, arcdist))
    return elev_angle

def get_bearing(lon1, lat1, lon2, lat2):
    """
    Gives the bearing clockwise from north in degrees from point 1 to
    point 2, at point 1 (init_bearing) and at point 2 (final_bearing)

    *** PARAMS ***
      lon1 [float]: longitude of point 1 (deg)
      lat1 [float]: latitude of point 1 (deg)
      lon2 [float]: longitude of point 2 (deg)
      lat2 [float]: latitude of point 2 (deg)

    ***RETURNS***
      init_bearing [float]: bearing clockwise from north in degrees at
        point 1's latitude and longitude.
      final_bearing [float]: final bearing after following ray along
        bearing from point 1 to point 2's latitude and longitude.
    """
    atan2_argx = np.sin(np.deg2rad(lon2 - lon1))*np.cos(np.deg2rad(lat2))
    atan2_argy1 = np.cos(np.deg2rad(lat1))*np.sin(np.deg2rad(lat2))
    atan2_argy2 = np.sin(np.deg2rad(lat1))*np.cos(np.deg2rad(lat2)) \
        * np.cos(np.deg2rad(lon2 - lon1))
    init_bearing = np.rad2deg(np.arctan2(atan2_argx, atan2_argy1-atan2_argy2))

    # The FINAL bearing is just the bearing from final point to the initial
    # point, reversed by adding 180 modulo 360
    atan2_argx = np.sin(np.deg2rad(lon1 - lon2))*np.cos(np.deg2rad(lat1))
    atan2_argy1 = np.cos(np.deg2rad(lat2))*np.sin(np.deg2rad(lat1))
    atan2_argy2 = np.sin(np.deg2rad(lat2))*np.cos(np.deg2rad(lat1)) \
        * np.cos(np.deg2rad(lon1 - lon2))
    rev_bearing = np.rad2deg(np.arctan2(atan2_argx, atan2_argy1-atan2_argy2))
    final_bearing = (rev_bearing + 180.) % 360.
    return init_bearing, final_bearing

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
    Return the rotation matrix associated with counterclockwise rotation
    about the given axis by theta radians.
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


def dir_ned2geo(loc_sph, dir_NED, time=None):
    """
    Convert a direction vector in terms of N-E-D components to the GEI
    coordinate system.

    *** PARAMS ***
        loc_sph: spherical geographic coordinates location:
            (alt, Altitude of point (km)
             lat, Latitude of point (deg)
             lon) Longitude of point (deg)
        dir_NED: pointing direction in NED coordinates
            (N, direction component in north direction
             E, direction component in east direction 
             D) direction component in down direction
        [time]: time for which this applies [default: None] [datetime]

    *** RETURNS ***
        dir_XYZ: direction component in cartesian GEO's X/Y/Z directions

    """
    from spacepy import coordinates as coord
    from spacepy import time as tm
    import datetime as dt
    if time is None:
        time = dt.datetime.now()

    (N, E, D) = (dir_NED[0], dir_NED[1], dir_NED[2])
    (alt, lat, lon) = (loc_sph[0], loc_sph[1], loc_sph[2])
    (lt_r, ln_r) = (np.deg2rad(lat), np.deg2rad(lon))
    dir_NED_arr = np.array(dir_NED)

    # ( c_x )   ( -cos(phi)sin(theta) -sin(theta) -cos(phi)cos(theta) )( c_n )
    # ( c_y ) = ( -sin(phi)sin(theta)  cos(theta) -sin(phi)cos(theta) )( c_e )
    # ( c_z )   (       cos(theta)         0          -sin(theta)     )( c_d )
    #
    # Where phi is longitude (or ln_r) and theta is latitude (lt_r)

    A_x = np.array(( -np.cos(ln_r)*np.sin(lt_r), -np.sin(lt_r), 
                     -np.cos(ln_r)*np.cos(lt_r) ))
    A_y = np.array(( -np.sin(ln_r)*np.sin(lt_r), np.cos(lt_r),
                     -np.sin(ln_r)*np.cos(lt_r) ))
    A_z = np.array(( np.cos(lt_r), 0, -np.sin(lt_r) ))
    A = np.array([A_x, A_y, A_z])
    dir_XYZ = np.dot(A,dir_NED_arr)
    #print("dir_XYZ: \n",dir_XYZ)
    
    logging.info("dir_XYZ: \n{0}".format(dir_XYZ))
    return dir_XYZ

# ----------------------------------------------------------------------------- 

def initialize_logger():
    """
    Function for setting up the initial logging parameters

    :param use_verbose: [boolean] flag indicating whether to be verbose.
        ** If _not_ running parse/fetch requests from the command-line **
    """
    global rLogger
    level = logging.INFO if quiet_mode else logging.DEBUG
    LOG_FILE = 'conjunctions.log'

    logging.basicConfig(level=level,
        format='%(levelname)s %(asctime)s: %(message)s', 
        datefmt='%m/%d/%Y %I:%M:%S %p')
 
    logFormatter = logging.Formatter('%(levelname)s %(asctime)s: %(message)s')
    rLogger = logging.getLogger(__name__)
    rLogger.setLevel(level)

    fileHandler = logging.FileHandler("./{0}".format(LOG_FILE))
    fileHandler.setFormatter(logFormatter)
    rLogger.addHandler(fileHandler)



# pylint: disable=C0103
if __name__=="__main__":
    """ TESTING INDEX OF REFRACTION-RELATED FUNCTIONS """
    freq10 = 1.0422E7   #[Hz]
    freq12 = 1.2500E7   #[Hz]
    txlon = -75.552    #[Deg]
    txlat = 45.403     #[Deg]

    date_string="20160421"
    filename, __ = data_utils.get_ottawa_data(date_string)
    omega_p = 2*np.pi*5.975E6 # [rad/s]
    omega_c = 2*np.pi*1.3E6   # [rad/s]
    fof2_alt = 260. # [km]

    lons, lats, alts, ephtimes = data_utils.get_rri_ephemeris(filename)
    b_igrf = magnet_data.get_igrf(lons, lats, alts, ephtimes)

    # Some extra stuff I want to do for analysis at the bottom here...
    # Load Rob Gillies' ephemeris he used for his ray trace plots
    Roblons, Roblats, Robalts = data_utils.load_rob_ephemeris(
                                data_utils.RRITK_DATA + '/satcoords_20160418.txt')
    # Load ephemeris from 20160418 for ray trace plots and plasma
    # intersection testing
    lons18, lats18, alts18, ephtimes18 = data_utils.get_rri_ephemeris(
        data_utils.get_ottawa_data('20160418')[0])

    # For validating get_plasma_intersection manually by plotting etc.
    # Validate get_plasma_intersection manually in OttawaPlots.ipynb
    plons=[]
    plats=[]
    i = 19
    lon18 = lons18[i]
    lat18 = lats18[i]
    alt18 = alts18[i]
    ephtime18 = ephtimes18[i]
    time18 = ephem_to_datetime(ephtime18)
    plalts = np.arange(60., alt18)
    for plalt in plalts:
        plon, plat = get_plasma_intersection(lon18, lat18, alt18, plasma_alt=plalt)
        plons.append(plon)
        plats.append(plat)
    
    # ALSO: Testing output of approximate aspect angle method vs rigorous:
    angs_apprx = []
    angs_rigor = [] 
    kvs = []
    kvs2 = []
    bvs = []
    bvs2 = []
    for i in range(len(alts18)):
        # Particular location point
        lon18 = lons18[i]
        lat18 = lats18[i]
        alt18 = alts18[i]
        ephtime18 = ephtimes18[i]
        time18 = ephem_to_datetime(ephtime18)
        # Get the different versions of the line-of-sight vector
        kv2 = get_kvec2(lon18, lat18, alt18)
        kvs2.append(kv2)
        kv = get_kvecs(lon18, lat18, alt18)
        kvs.append(kv)
        # Get the different versions of the magnetic field vector
        bv2 = get_bvec2(lon18, lat18, alt18, time18)
        bvs2.append(bv2)
        bv = get_bvec(lon18, lat18, alt18, time18)
        bvs.append(bv)
        # Get the different ve
        ang_rigor = get_aspect_angle(kv2, bv2)
        angs_rigor.append(ang_rigor)
        ang_apprx = get_aspect_angle(kv, bv)
        angs_apprx.append(ang_apprx)
    
    datarr, datlats = data_utils.load_density_profile(
                        data_utils.RRITK_DATA + '/20160418-densities.txt')  
    tecs, mean_bcs, ql_integrals, phase_integrals = faraday_pass(
                lons18, lats18, alts18, ephtimes18, datarr, datlats)
    plt.plot(phase_integrals); plt.show()
    
