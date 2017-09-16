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
import datetime as dt

import numpy as np 
import matplotlib.pyplot as plt
import logging

from davitpy.utils import plotUtils

import __init__
import rriephtk.utils.data_utils as datautils
from rriephtk.utils.data_utils import two_pad as tp

OTTAWA_TX_LON = -75.552
OTTAWA_TX_LAT = 45.403
OTTAWA_TX_ELEV = 0.070 # 70m elevation in Ottawa - small compared to 300km satellite
MILLSTONE_TX_LON = -71.491 # acquired from google maps satellite imagery
MILLSTONE_TX_LAT = 42.619
data_path = '../rri-conjunction-script/data'

def plot_sat_ephemeris(date_string=None,lons=None,lats=None,alts=None,ephtimes=None, tx_lon=None, tx_lat=None):
    """
    Attempt to give a nice ground-track 

    """
    if lons is None or lats is None or alts is None or ephtimes is None:
        if date_string==None or not isinstance(date_string, type("e.g.")):
            print "Need to provide either lons/lats/alts/ephtimes or date_string!"
            return -1
        # e.g. like paired with the MGF data access
        rri_fname,indx_rev = datautils.get_ottawa_data(date_string) 
        lons,lats,alts,ephtimes = datautils.get_rri_ephemeris(rri_fname)
    times=  datautils.ephems_to_datetime(ephtimes)

    datautils.update_progress(0.1)    

    # A different font for the legend etc. might be nice
    font = {'fontname':'Computer Modern'}
    m = plotUtils.mapObj(lat_0=np.mean(lats), lon_0=np.mean(lons), width=4.0*(max(lons) - min(lons))*1000*180, \
                         height=1.3*(max(lats) - min(lats))*1000*180, coords='geo',resolution='i', \
                         datetime=times[0])
    # (the 1000* factors are to replace 1000 in how usually these are written as "width=111e3*180")

    datautils.update_progress(0.3)

    x,y = m(lons, lats, coords='geo')
    m.plot(x, y, 'b-', label="ePOP ground track")

    if tx_lon is not None and tx_lat is not None:
        x, y = m(tx_lon, tx_lat, coords='geo')
        m.plot(x, y, 'ro', label='Transmitter')        

    datautils.update_progress(0.6)

    plt.xlabel('Geographic Longitude (degrees)')
    plt.ylabel('Geographic Latitude (degrees)')
    times_str = tp(times[0].hour) + ":" + tp(times[0].minute) + \
                 "-" + tp(times[-1].hour) + ":" + tp(times[-1].minute)
    dates_str = str(times[0].year) + tp(times[0].month) + tp(times[0].day)
    plt.title("CASSIOPE Ground Track for " + dates_str + " " + times_str)
    #plt.legend(loc='best')
    datautils.update_progress(0.9)
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

def plot_fan_rri(lons, lats, alts, time, codes=['sas'], save=False):
    """
    Plot the RRI ground track along with specified SuperDARN power plots
    from the corresponding time period.
    """
    from davitpy import pydarn
    import fanrri
    assert(type(time) == dt.datetime)

    for code in codes:
        logging.info("Processing station '{0}'...".format(code)) 

    if save:
        fanrri.plotFan(time, codes, param='power', gsct=False, lons=lons, lats=lats, 
                        png=True, pdf=False, show=False) 
    else:
        fanrri.plotFan(time, codes, param='power', gsct=False, lons=lons, lats=lats, 
                        png=False, pdf=False, show=True) 
    #pydarn.plotting.fan.plotFan(time, codes, param='power', gsct=False)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
#               Plots for Danskin's 2017 Paper, Figure 1
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -

def plot_all5(geographic=False, latspacing=10.):
    """
    Plot April 18th-22nd ephemeris tracks on same mapobj
    """
    from rriephtk.utils.data_utils import *
    from rriephtk.analysis.analysis_tools import * 
    # Just in case we wanted to show points where we thought the ellipticity
    # angle experienced a distinctive shift also, these 'ellip_rev_xth' variables
    # are meant to show the index which we figure corresponds to their distinct 
    # shift (checking to see if it had any relationship to the faraday rotation
    # inflection point) 
    fname_18th, idx_rev_18th = get_ottawa_data("20160418")
    lons_18th, lats_18th, alts_18th, ephtimes_18th = get_rri_ephemeris(fname_18th)
    fname_19th, idx_rev_19th = get_ottawa_data("20160419")
    lons_19th, lats_19th, alts_19th, ephtimes_19th = get_rri_ephemeris(fname_19th)
    fname_20th, idx_rev_20th = get_ottawa_data("20160420")
    lons_20th, lats_20th, alts_20th, ephtimes_20th = get_rri_ephemeris(fname_20th)
    fname_21st, idx_rev_21st = get_ottawa_data("20160421")
    lons_21st, lats_21st, alts_21st, ephtimes_21st = get_rri_ephemeris(fname_21st)
    fname_22nd, idx_rev_22nd = get_ottawa_data("20160422")
    lons_22nd, lats_22nd, alts_22nd, ephtimes_22nd = get_rri_ephemeris(fname_22nd)
    
    times_18th = ephems_to_datetime(ephtimes_18th)
    times_19th = ephems_to_datetime(ephtimes_19th)
    times_20th = ephems_to_datetime(ephtimes_20th)
    times_21st = ephems_to_datetime(ephtimes_21st)
    times_22nd = ephems_to_datetime(ephtimes_22nd)

    indx_shortest_18th, dists_18th = get_closest_approach(lons_18th, lats_18th, alts_18th)
    indx_shortest_19th, dists_19th = get_closest_approach(lons_19th, lats_19th, alts_19th)
    indx_shortest_20th, dists_20th = get_closest_approach(lons_20th, lats_20th, alts_20th)
    indx_shortest_21st, dists_21st = get_closest_approach(lons_21st, lats_21st, alts_21st)
    indx_shortest_22nd, dists_22nd = get_closest_approach(lons_22nd, lats_22nd, alts_22nd)
    
    # The numeric data type that I was retrieving from geog_longs, when _NOT_ stored
    # in an array, was being rejected by the mapObj() function below. So I convert 
    # these numbers to floats explicitly here.
    shlon_18th = float(lons_18th[indx_shortest_18th])
    shlat_18th = float(lats_18th[indx_shortest_18th])
    invlon_18th = float(lons_18th[idx_rev_18th])
    invlat_18th = float(lats_18th[idx_rev_18th])
 
    shlon_19th = float(lons_19th[indx_shortest_19th])
    shlat_19th = float(lats_19th[indx_shortest_19th])
    invlon_19th = float(lons_19th[idx_rev_19th])
    invlat_19th = float(lats_19th[idx_rev_19th])

    shlon_20th = float(lons_20th[indx_shortest_20th])
    shlat_20th = float(lats_20th[indx_shortest_20th])
    invlon_20th = float(lons_20th[idx_rev_20th])
    invlat_20th = float(lats_20th[idx_rev_20th])

    shlon_21st = float(lons_21st[indx_shortest_21st])
    shlat_21st = float(lats_21st[indx_shortest_21st])
    invlon_21st = float(lons_21st[idx_rev_21st])
    invlat_21st = float(lats_21st[idx_rev_21st])
 
    shlon_22nd = float(lons_22nd[indx_shortest_22nd])
    shlat_22nd = float(lats_22nd[indx_shortest_22nd])
    invlon_22nd = float(lons_22nd[idx_rev_22nd])
    invlat_22nd = float(lats_22nd[idx_rev_22nd])
   
    # A different font for the legend etc. might be nice
    fig = plt.figure(1, figsize=(14,8))
    font = {'fontname':'Computer Modern'}

    if geographic==True:
        m = plotUtils.mapObj(lat_0=46.0, lon_0=-75.0, width=18e3*180, height=25e3*90, coords='geo',resolution='i',datetime=times_20th[0], gridLatRes=latspacing, fillOceans=(1.,1.,1))
    else:
        m = plotUtils.mapObj(lat_0=57.0, lon_0=5.0, width=18e3*180, height=25e3*90, coords='mag',resolution='i',datetime=times_20th[0], gridLatRes=latspacing, fillOceans=(1.,1.,1))

    # FIRST: Plot the location of Ottawa
    x,y = m(OTTAWA_TX_LON,OTTAWA_TX_LAT,coords='geo')
    m.plot(x,y,'r-o',markersize=8,label="Ottawa")
    x,y = m(MILLSTONE_TX_LON,MILLSTONE_TX_LAT,coords='geo')
    m.plot(x,y,'m-o',markersize=8,label="Millstone Hill Digisonde")
    
    # SECOND: Plot the satellite ground-track.
    x,y = m(lons_18th, lats_18th, coords='geo')
    m.plot(x,y,'b-',label="18 April")#label="EPOP ground track")
    x,y = m(lons_19th, lats_19th, coords='geo')
    m.plot(x,y,'b:',label="19 April")
    x,y = m(lons_20th, lats_20th, coords='geo')
    m.plot(x,y,'b--',label="20 April")
    x,y = m(lons_21st, lats_21st, coords='geo')
    m.plot(x,y,'k-',label="21 April")
    x,y = m(lons_22nd, lats_22nd, coords='geo')
    m.plot(x,y,'k--',label="22 April")
    
    # THIRD: Plot a circle emphasizing the point of closest approach
    x,y = m(shlon_18th, shlat_18th, coords='geo')
    m.plot(x,y,'bo',markersize=6,label='Closest Approach TX')
    x,y = m(shlon_19th, shlat_19th, coords='geo')
    m.plot(x,y,'bo',markersize=6)
    x,y = m(shlon_20th, shlat_20th, coords='geo')
    m.plot(x,y,'bo',markersize=6)
    x,y = m(shlon_21st, shlat_21st, coords='geo')
    m.plot(x,y,'bo',markersize=6)
    x,y = m(shlon_22nd, shlat_22nd, coords='geo')
    m.plot(x,y,'bo',markersize=6)

    # FOURTH: Plot the point I've determined is the point of the Faraday Rotation inversion.
    x,y = m(invlon_18th, invlat_18th, coords='geo')
    m.plot(x,y,'go',markersize=6,label=(r'Reversal Point of Faraday Rotation $(\psi)$'))
    x,y = m(invlon_19th, invlat_19th, coords='geo')
    m.plot(x,y,'go',markersize=6)
    x,y = m(invlon_20th, invlat_20th, coords='geo')
    m.plot(x,y,'go',markersize=6)
    x,y = m(invlon_21st, invlat_21st, coords='geo')
    m.plot(x,y,'go',markersize=6)
    x,y = m(invlon_22nd, invlat_22nd, coords='geo')
    m.plot(x,y,'go',markersize=6)

    ax = plt.gca()
    #plt.xlabel('Magnetic Longitude (degrees)')
    if geographic==True:
        ax.set_xlabel('Geographic Longitude (degrees)')
        plt.ylabel('Geographic Latitude (degrees)')
    else: 
        ax.set_xlabel('Magnetic Longitude (degrees)')
        plt.ylabel('Magnetic Latitude (degrees)')
    ax.xaxis.set_label_coords(0.5,-0.050)
    
    #plt.title("ePOP Pass vs. Ottawa radar 18-22 April 2016")
    plt.legend(loc='lower right', numpoints = 1,fontsize=11)#loc='best')
    plt.savefig('tmp_all5.eps', format='eps', bbox_inches='tight')
    plt.savefig('tmp_all5.png', format='png', bbox_inches='tight')
    plt.show()

def do_fig1_plots():
    """
    There were a few different options we tried before settling on how 
    figure 1 should be    
    """
    plot_all5(geographic=True, latspacing=10.)
    plot_all5(geographic=True, latspacing=5.)   
    plot_all5(geographic=False, latspacing=10.)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -


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
