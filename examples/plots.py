#!/usr/bin/env python
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

import matplotlib
import mpl_toolkits
import mpl_toolkits.basemap

import numpy as np 
import matplotlib.pyplot as plt
import argparse

from davitpy.utils import plotUtils

import rriephtk_path
import rriephtk.utils.data_utils as datautils
from rriephtk.utils.data_utils import two_pad as tp
from rriephtk import plotting

#import logging
#logging.basicConfig(level=logging.DEBUG,
#    format='%(levelname)s %(asctime)s: %(message)s', 
#    datefmt='%m/%d/%Y %I:%M:%S %p')
#logger = logging.getLogger(__name__)


 
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -

def get_args():
    """
    Parse the command-line arguments.

    Yes, in an ideal world, this whole thing would be a sweet little 
    object which does things on initialization, but at least for now,
    this works as a stand-alone function!

    ** RETURNS **
    rrifile (string): filename
    rcode (string): radar code as 3-letter string (e.g. 'sas')
    save_plot (string): user directive on whether or not to output plot as a saved .png
    """
    parser = argparse.ArgumentParser()

    # For now, we require a particular station to be requested
    parser.add_argument("-r", "--radar_code", 
                        help="SuperDARN Station code to produce fan plot for (e.g. sas)",
                        type=str)

    parser.add_argument("-f", "--rri_file", help="Specified RRI file to plot ephemeris",
                        type=str)

    parser.add_argument("-s", "--save_plot", help="Whether to save the plot instead of showing (y/n)",
                        action="store_true", default=False)
    args = parser.parse_args()
    rrifile = args.rri_file
    rcode = args.radar_code
    save_plot = args.save_plot
    return rrifile, rcode, save_plot

def process_args(rrifile, rcode, save_plot):
    """
    Encapsulates the necessary logic to decide what to do based on 
    command-line arguments.
    """
    if rrifile is not None:
        if rcode is not None:
            print("Attempting to use 'plot_fan_rri'")
            #TODO: validate input
            try:
                plotting.plot_fan_rri(lons, lats, alts, ephtimes, codes=[rcode], save=save_plot)
            except Exception as e:
                print("Problem plotting fanplot: {0}".format(e))            
        else:
            print("Attempting to print the satellite ephemeris nicely")
            try:
                lons, lats, alts, ephtimes = datautils.get_rri_ephemeris(rrifile)
            except Exception as e:
                print("Problem opening the RRI file - aborting! Exception was: {0}".format(e))
                return None
            params = plot_options()
            plotting.plot_sat_ephemeris(lons, lats, alts, ephtimes, save=save_plot, 
                w_mult=params["w_mult"], h_mult=params["h_mult"], lon_offs=params["lon_offs"], 
                lat_offs=params["lat_offs"])

    return None 

def plot_options():
    """
    Function to do the queries for additional information for plot_sat_ephemeris

    :params: None
    :return: [dict] options: dictionary of values for w_mult (scalar mult), h_mult (scalar mult), 
                            long_offs (deg), lat_offs (deg)
    """
    offs = dict()
    mults = dict()
    options = dict()
    print("Multiplicative factor for window width [Default 1.0]: ")
    mults["w_mult"] = raw_input()
    print("Multiplicative factor for window height [Default 1.0]: ")
    mults["h_mult"] = raw_input()
    print("Longitude offset for plot center with respect to sat track [Default 0 deg]: ")
    offs["lon_offs"] = raw_input()
    print("Latitude offset for plot center with respect to sat track [Default 0 deg]: ")
    offs["lat_offs"] = raw_input()
    
    for k in offs.keys():
        options[k] = offs[k] if (type(offs[k]) in [int, float] and offs[k] >= 0) else 0.
    for k in mults.keys():
        options[k] = mults[k] if (type(mults[k]) in [int, float] and mults[k] > 0) else 1. 
    return options
    

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -

if __name__=="__main__":
    fname, rcode, save_plot = get_args()
    #path, fname = datautils.initialize_data()
    #lons, lats, alts, ephtimes = datautils.get_rri_ephemeris(fname)
    #plot_fov_sat("Saskatoon", lons, lats, date=datautils.ephem_to_datetime(ephtimes[0])) 
    print("Args were\tfname: {0},\trcode: {1},\tsave_plot: {2}".format(fname, rcode, save_plot))
    process_args(fname, rcode, save_plot)
