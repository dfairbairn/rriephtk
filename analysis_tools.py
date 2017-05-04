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

import numpy as np
import matplotlib.pyplot as plt
import data_utils
import magnet_data
import ottawa_plots
""" 
Index of Refraction-related Code
--------------------------------
    The purpose of this code is the calculation of index of refraction vs
    incidence angle for radio wave propagation direction compared to the
    ionospheric plasma's orientation/local B field

"""

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

def cyclotron_freq(B):
    """
    Takes a magnetic field value(s) parameter and uses it to calculate the
    cyclotron frequency for an electron in this field.

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
        angles [float]: Aspect angles (computed by get_kb_angle()
        nplus [float]: positive Appleton-Hartree solution (O-mode index of ref)
        nminus [float]: negative AH solution (X-mode index of ref)

    """
    from ottawa_plots import get_kb_ottawa_angle
    # TODO: generalize get_kb_ottawa_angle to get_kb_angle
    # kvecs,bvecs,angles = get_kb_angle(lons,lats,alts,ephtimes,tx_lon,tx_lat)
    bvecs,kvecs,angles = get_kb_ottawa_angle(lons,lats,alts,ephtimes)
    angles_rad = np.deg2rad(angles)
    #TODO: allow user to input nu_e to this?
    X,Y,Z = appleton_coeffs(omega_c,omega_p,freq,0.0) 
    np_sq,nm_sq = appleton_hartree(X,Y,Z,angles_rad)
    nplus = np.sqrt(np_sq)
    nminus = np.sqrt(nm_sq)
    return angles_rad, nplus, nminus 

if __name__=="__main__":

    """ TESTING INDEX OF REFRACTION-RELATED FUNCTIONS """
    freq10 = 1.0422E7   #[Hz]             
    freq12 = 1.2500E7   #[Hz]
    tx_lon = -75.552    #[Deg]
    tx_lat = 45.403     #[Deg]
 
    filename, __ = ottawa_plots.get_ottawa_data("20160421")
    omega_p = 2*np.pi*5.975E6 # [rad/s]
    fof2_alt = 260. # [km] 
    
    lons,lats,alts,ephtimes = data_utils.get_rri_ephemeris(filename)
    b_igrf = magnet_data.get_igrf(lons,lats,alts,ephtimes)

    # Test improved cyclotron frequency function
    omega_c1 = improved_cyclotron_freq(lons,lats,alts,ephtimes,fof2_alt=280.)

    # Test regular cyclotron frequency function
    omega_c2 = cyclotron_freq(b_igrf)

    # Test basic_ionosphere_params
    basic_omega_c, basic_omega_p, basic_l_d = basic_ionosphere_params(altitude=300.)

    # Test appleton_coeffs
    X,Y,Z = appleton_coeffs(omega_c1, omega_p, freq10, 0.)    

    # Test appleton_hartree
    theta = np.array( range(315) )/100.
    nplus,nminus = appleton_hartree(X,Y,Z,theta)

    # Test txpass_to_indices 
    angles_r,nplus,nminus = txpass_to_indices(lons,lats,alts,ephtimes,tx_lon,tx_lat,freq10,omega_p,omega_c1)
