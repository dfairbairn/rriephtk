# -*- coding: utf-8 -*-
"""
file: 'test_antools.py'
description:
    This document contains the unit tests for analysis_tools.
author: David Fairbairn
date: May 2017

"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import logging

import rriephtk.utils.data_utils as data_utils
import rriephtk.utils.magnet_data as magnet_data

from rriephtk.utils.data_utils import ephem_to_datetime, ephems_to_datetime 
from rriephtk.analysis.analysis_tools import *

OTTAWA_TX_LON = -75.552
OTTAWA_TX_LAT = 45.403
OTTAWA_TX_ELEV = 0.070 # Ottawa 70m elevation - small compared to satellite

# Physical constants:
EL_CHARGE = 1.602E-19   #[C]
EL_MASS = 9.109E-31 #[kg]
EPS0 = 8.8542E-12   #[A*s/(V*m)]
EARTH_RAD = 6371.   #[km]

logging.basicConfig(filename='./data/analysis-tools.log',level=logging.DEBUG)

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

    '''
    I'm primarily checking that the functions don't break in general, but the
    cases below check particular numbers in order to detect whether the data or
    the processing unexpectedly changes after some updates to things.
    '''
    omega_c1 = improved_cyclotron_freq(lons, lats, alts, ephtimes, fof2_alt=280)
    omega_c2 = cyclotron_freq(b_igrf)
    print("Beginning unit testing; be forewarned," \
        + " IGRF calls will plug up standard out...")

    # Test regular and improved cyclotron frequency functions
    # ------------------
    if np.round(np.log10(omega_c2)) != 7 or np.round(np.log10(omega_c1)) != 7:
        logging.error("Error with cyclotron frequency calculations")

    # Test plasma_freq(n_e)
    # ------------------
    pf = plasma_freq(5E11)
    if np.round(np.log10(pf),3) != 7.601:
        logging.error("Error with plasma_freq()")

    # Test basic_ionosphere_params
    # ------------------
    w_c, w_p, l_d = basic_ionosphere_params(altitude=300.)
    if np.round(np.log10(w_p),1)!=7.8 or np.round(np.log10(w_c),1) != 7.0:
        logging.error("Error with basic ionosphere params")

    # Test appleton_coeffs
    # ------------------
    # X should be about 0.25-0.3, Y about 0.1, Z about 0
    X, Y, Z = appleton_coeffs(omega_c1, omega_p, freq10, 0.)
    if np.round(X,2) != 0.33 or np.round(Y,2) != 0.13 or np.round(Z) != 0.:
        logging.error("Error with appleton_coeffs()")

    # Test appleton_hartree
    # ------------------
    theta = np.array(range(315))/100.
    npl1, nmin1 = appleton_hartree(X, Y, Z, theta)
    if np.round(abs(npl1[0]),3)!=0.842 or np.round(abs(nmin1[0]),3)!=0.790:
        logging.error("Error with appleton_hartree()")

    '''
    # Test txpass_to_indices
    # ------------------
    angs_r, npl2, nmin2 = txpass_to_indices(lons, lats, alts, ephtimes, txlon,
                                            txlat, freq10, omega_p, omega_c1)
    if np.round(np.mean(angs_r),3)!=2.320 or np.round(abs(npl2[0]),3)!=0.914:
        logging.error("Error with txpass_to_indices()")
    '''

    """ TESTING DIRECTIONAL/EPHEMERIS-RELATED FUNCTIONS"""

    # Test get_bearing()
    # ------------------
    # Correctness:
    tst_bearing1,__ = get_bearing(-106.,52.,-96.,52.)
    tst_bearing2,__ = get_bearing(-106.,52.,-106.,42.)
    if np.round(tst_bearing1,1) != 86.1 or np.round(tst_bearing2,1) != 180.:
        logging.error("Error calculating bearings")
    # Parallelization:
    a = np.array([100.,100.,100.])
    b = np.array([45.,45.,45.])
    c = np.array([35.,55.,-45.])
    bearings_i, bearings_f = get_bearing(a, b, a, c)
    if (bearings_i != [180.,0.,180.]).any(): # if any entries incorrect...
        logging.error("Error with parallelized bearings calculation")

    # Test haversine:
    # ------------------
    hav = haversine(50,50,60,60)
    if int(hav)!=1278: #corroborated with other calculators
        logging.error("Error with haversine()")
    # check parallelization
    havs = haversine(a, b, a, c)

    # Test elevation_angle:
    # ------------------
    elev_a = get_elevation_angle(-106.,52.,0.,-96.,52.,300.)
    if np.round(elev_a,2) != 23.68:
        logging.error("Error with get_elevation_angle()")

    # Test rotation_matrix:
    # ------------------
    rotA = np.round(rotation_matrix([1,0,0],-np.pi/6.),3)
    rotB = np.round(rotation_matrix([1,0,0],2*np.pi),3)
    answerA = np.round(np.array([[1.,0.,0.],[0., np.sqrt(3.)/2.,1/2.],
                       [0,-1/2., np.sqrt(3.)/2.]]),3)
    answerB = np.round(np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]),3)
    if (rotA != answerA).any() or (rotB != answerB).any():
        logging.error("Error with rotation_matrix()")

    # Test get_kvecs:
    # ------------------
    kvs = get_kvecs(lons, lats, alts) # From Ottawa by default.
    # I've already checked the validity of the April 21st data.
    # This test relies on it:
    if (np.round((np.mean(kvs[:,0]), np.mean(kvs[:,1]), np.mean(kvs[:,2])),2) \
            != [0.09,0.05,-0.65]).any():
        logging.error("Error with get_kvecs()")

    # Test get_igrf(): already called get_igrf earlier in this if__name__ block
    # ------------------
    #b_igrf = get_igrf(lons, lats, alts, ephtimes)
    if (np.round(b_igrf[0]) != np.array([17858.,-3148.,39417.])).any():
        logging.error("Error with get_igrf()")

    # Test get_aspect_angle():
    # ------------------
    a1 = np.array((1, 0, 0))
    b1 = np.array((1, 1, 0))
    c1 = np.array((-1, 0, 1))
    kv_tst = np.array((1,2,3))
    bv_tst = np.array((2,30,7))
    tst1 = np.round(get_aspect_angle(a1, b1), 1)==45.0
    tst2 = np.round(get_aspect_angle(a1, c1), 1)==135.0
    tst3 = np.round(get_aspect_angle(kv_tst,bv_tst), 2)==44.06
    if not tst1 or not tst2 or not tst3:
        logging.error("Error with get_aspect_angle()")

    # Test get_plasma_intersection():
    # ------------------
    # Use a point directly north, check that intersections are also north
    plasma_lon, plasma_lat = get_plasma_intersection(-75.552,48.403,370.)
    if np.round(plasma_lon,1)!=-75.6 or np.round(plasma_lat,1)!=47.8:
        logging.error("Error with get_plasma_intersection")

    # Test get_kb_angle():
    # ------------------
    # Already tested successfully earlier by doing tx_pass_indices
    #bvs, kvs, angles_r = get_kb_angle

    # Test get_kdip_angles():
    # ------------------
    #TODO: test and *validate*!
    (lons, lats, alts, ephtimes, mlons, mlats, mlts, pitch, yaw, roll) = \
        data_utils.get_rri_ephemeris_full(filename)
    dipole_dirs, kdip_angles = get_kdip_angles(lons, lats, alts, ephtimes,
                                               pitch, yaw, roll)

    # Test get_ramdirs():
    # ------------------
    vs, dists = get_ramdirs(lons, lats, alts)
    if (np.round(vs[0][0],3)!=7.442) or (np.round(dists[-1],3)!=7.419):
        logging.error("Error with get_ramdirs")

    # Test get_closest_approach:
    # ------------------
    index_closest, dists = get_closest_approach(lons, lats, alts)
    if (index_closest!=103) or np.round(dists[-1],2)!=1049.97:
        logging.error("Error with get_closest_approach")
    """
    # Test faraday_trace
    # ------------------
    datarr, datlats = data_utils.load_density_profile(
                        './data/20160418-densities.txt')
    tec, bcs, ql_int, faraday_int = faraday_trace(lons[0], lats[0], alts[0],
                                               ephtimes[0], datarr, datlats)
    # validated calculations by noting similarity of numbers with Rob's
    if np.round(np.log10(tec),3)!=17.198:
        logging.error("Error with faraday_trace")
    """

    # Some extra stuff I want to do for analysis at the bottom here...
    # Load Rob Gillies' ephemeris he used for his ray trace plots
    Roblons, Roblats, Robalts = data_utils.load_rob_ephemeris(
                                './data/satcoords_20160418.txt')
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

  
    # Test dir_ned2geo 
    # ----------------
    from spacepy import coordinates as coord
    from spacepy import time as tm
    import datetime as dt
    
    # Test along the 0 degrees meridian (varies in X, Z directions)
    aas = []
    bbs = []
    for l in np.arange(-90., 90.):
        a = coord.Coords([0., l, 0.],'GEO','sph')
        aas.append(a)
        a.ticks = tm.Ticktock(dt.datetime.now())
        b = a.convert('GEO', 'car')
        bbs.append(b)

    # Radially inward vector at each prime meridian pt (N, E, D) = (0, 0, 1)
    # Should go (0, 0, 1) to (-1, 0, 0) at 0 latitude to (0, 0, -1)
    dirs_geo = [ dir_ned2geo(aa.data[0], (0, 0, 1)) for aa in aas ]
    tst1 = (np.round(dirs_geo[0], 1)==(0., 0., 1.)).all()
    tst2 = (np.round(dirs_geo[len(dirs_geo)/2], 1)==(-1., 0., 0.)).all()
    tst3 = (np.round(dirs_geo[-1], 1)==(0., 0., -1.)).all()
    if not tst1 or not tst2 or not tst3:
        logging.error("Error with dirs_ned2geo()")

    # Test along the 90 degrees meridian (varies in Y, Z directions)
    ccs = []
    dds = []
    for l in np.arange(-90., 90.):
        c = coord.Coords([0., l, 90.],'GEO','sph')
        ccs.append(c)
        c.ticks = tm.Ticktock(dt.datetime.now())
        d = a.convert('GEO', 'car')
        dds.append(d)
    # 'Northward' facing vector at each 90 deg meridian pt (N, E, D) = (1, 0, 0)
    # Should go (0, 1, 0) to (0, 0, 1) to (0, -1, 0)
    dirs_geo_2 = [ dir_ned2geo(cc.data[0], (1, 0, 0)) for cc in ccs ]
    tst1 = (np.round(dirs_geo_2[0], 1)==(0., 1., 0.)).all()
    tst2 = (np.round(dirs_geo_2[len(dirs_geo_2)/2], 1)==(0., 0., 1.)).all()
    tst3 = (np.round(dirs_geo_2[-1], 1)==(0., -1., 0.)).all()
    if not tst1 or not tst2 or not tst3:
        logging.error("Error with dirs_ned2geo() bogooibiobjboj")

    # Test get_kvec2
    # ---------------
    # Take the coords along the prime meridian 'p' and just take the triplets
    p = [aa.data[0] for aa in aas]
    # Convert point-to-point directions (each is 'north')
    # get_kvec2(<lon>, <lat>, <alt>, <tx_lon= >, <tx_lat= >, <tx_alt= >)
    kvs = [ get_kvec2(p[i][2], p[i][1], p[i][0], tx_lon=p[i-1][2],
            tx_lat=p[i-1][1], tx_alt=p[i-1][0]) for i in np.arange(1, len(p))]
    # Expect a northward dir at south pole =+x dir, at eq =+z, at north pole =-x
    tst1 = (np.round(kvs[0],1)==(1.,0.,0.)).all()
    tst2 = (np.round(kvs[-1],1)==(-1.,0.,0.)).all()
    tst3 = (np.round(kvs[len(kvs)/2],1)==(0.,0.,1.)).all()
    if not tst1 or not tst2 or not tst3:
        logging.error("Error with get_kvec2()")

    # Test get_bvec2
    # ---------------
    # Hard to test this function, as it basically just takes the get_bvec
    # output and then converts it.

    # For now I'll try comparing the angle between get_bvec2's vectors
    # and get_kvec's vectors *after* they've been fed through the
    # dir_ned2geo function. These should be the same as 'angs' below.
    
    # I have independently confirmed that bv and bv2 are equal below:
    # bv = get_bvec(lon18, lat18, alt18, time18)
    # bv_conv = dir_ned2geo((alt18, lat18, lon18), bv)
    # bv2 = get_bvec2(lon18, lat18, alt18, time18)

    """    
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
    """ 
    """ 
    tecs, mean_bcs, ql_integrals, phase_integrals = faraday_pass(
                lons18, lats18, alts18, ephtimes18, datarr, datlats)
    plt.plot(phase_integrals); plt.show()
    """


    print("Tests complete!")
