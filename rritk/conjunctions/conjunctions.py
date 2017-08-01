"""
file: 'conjunctions.py'
description:
    This is the 'script 2.0' file. The functions herein allow 
    determination of intersections in CASSIOPE ephemeris during which
    it passes through the field-of-view of a SuperDARN radar.

    Basic functions simply report the intersection, and additionally
    the .h5 data files can have ephemeris points tagged with the 
    intersection details, and sections from the SuperDARN data logs can
    even be retrieved.

author: David Fairbairn
date: May 2017

This is a reworking/improvement of the original 'script.py' from 2016.
"""

import os
import subprocess
import sys
import logging

import davitpy
from davitpy import pydarn
import timeit
import math

import datetime as dt
import numpy as np

import __init__
import rritk.utils.data_utils as data_utils
import rritk.utils.range_cells as range_cells
from rritk.utils.data_utils import two_pad

OUTPUT_DIR = data_utils.RRITK_OUTPUT

class RRISuperdarnConjunction:
    """
    An object for representing a SuperDARN and ePOP RRI conjunction,
    AKA a case in which the CASSIOPE satellite's RRI instrument is
    recording measurements while also in proximity to SuperDARN radars.
    """

    def __init__(self, code, name, forb, bm_start, gt_start, bm_end, gt_end,
                 idx_start=-1, idx_end=-1):
        self.code = code
        self.name = name
        self.forb = forb
        self.bm_start = bm_start
        self.bm_end = bm_end
        self.gt_start = gt_start
        self.gt_end = gt_end
        self.idx_start = idx_start
        self.idx_end = idx_end
    
    def __repr__(self):
        str1 = "RRISuperDARNConjunction: '{0}' ({1}) ".format(self.code, self.forb)
        str2 = "Beam,Gate {0} to {1}".format((self.bm_start, self.gt_start), 
                                              (self.bm_end, self.gt_end))
        return str1+str2

    def conj_heuristic(self):
        """
        Judges the value of the conjunction in terms of how strongly RRI 
        could be affected by the SuperDARN radar.

        Proximity and extent of FOV crossed are the considerations.

        Note: 1000km proximity is 22 range gates of less. 
              1800km = 40 range gates
        """
        conj_extent = abs(self.bm_end - self.bm_start)
        conj_closeness = min(self.gt_end, self.gt_start)
        # Magic formula: '2' for large extent that's close. '1' for large 
        # extent or close, '0' for distant and small extents.
        score = 0
        if conj_extent >= 2.0:
            score += 1
        if conj_closeness <= 40.0:
            score += 1
        return score


def get_conjunctions(fname):
    """
    Returns a list of conjunction objects describing how the CASSIOPE
    ground-track intersects the field of view of SuperDARN radars. 

    """
    lons, lats, alts, ephtimes = data_utils.get_rri_ephemeris(fname)
    times = data_utils.ephems_to_datetime(ephtimes)

    nw = pydarn.radar.network()
    results = dict()
    # Subsets of relevant things
    lat_subs = lats
    lon_subs = lons

    relevant_radars = dict()
    conjs = []

    start_t = timeit.default_timer()
    activerads = [rad for rad in nw.radars if rad.status==1]

    for i, rad in enumerate(activerads):
        # Check active radars
        if rad.status == 1:
            fov_f = pydarn.radar.radFov.fov(site=rad.sites[0], altitude=300.,
                model='IS', coords='geo', ngates=75, fov_dir='front')
            bm_f, gt_f = range_cells.findRangeCell(lat_subs, lon_subs, fov_f)

            fov_b = pydarn.radar.radFov.fov(site=rad.sites[0], altitude=300.,
                model='IS', coords='geo', ngates=75, fov_dir='back')
            bm_b, gt_b = range_cells.findRangeCell(lat_subs, lon_subs, fov_b)
            
            results[rad.name] = (bm_f, gt_f, bm_b, gt_b)
            # non_nan_X contains indices of items that aren't nan (ie. within fov)
            non_nan_f = [n for n in range(np.size(bm_f)) if not math.isnan(bm_f[n])]        
            non_nan_b = [n for n in range(np.size(bm_b)) if not math.isnan(bm_b[n])]
             
            if (non_nan_b.__len__() > 0):
                rLogger.debug("Back field-of-view conjunction for {0}".format(rad))
                start = non_nan_b[0] 
                end = non_nan_b[ non_nan_b.__len__() - 1 ]
                conjs.append(RRISuperdarnConjunction(
                    rad.code[0], rad.name, 'back', bm_b[start], gt_b[start], 
                    bm_b[end], gt_b[end], idx_start=start, idx_end=end))
                relevant_radars[rad.name] = (rad.code[0],"back", bm_b[start],gt_b[start],bm_b[end],gt_b[end])

            if (non_nan_f.__len__() > 0):
                rLogger.debug("Front field-of-view conjunction for {0}".format(rad))
                start = non_nan_f[0]
                end = non_nan_f[ non_nan_f.__len__() - 1 ]
                conjs.append(RRISuperdarnConjunction(
                    rad.code[0], rad.name, 'front',  bm_f[start], gt_f[start], 
                    bm_f[end], gt_f[end], idx_start=start, idx_end=end))
                relevant_radars[rad.name] = (rad.code[0],"front",bm_f[start],gt_f[start],bm_f[end],gt_f[end])
        
        #Show progress on screen
        sys.stdout.flush()
        data_utils.update_progress((i+1)/float(len(activerads)))

    end_t = timeit.default_timer()
    rLogger.info("\nTime req'd to compute detailed intersections by brute force: " + \
        str(end_t - start_t) + " seconds.")

    # In general there will be a post-processing step in which only the 
    # conjunction results will be included in the output
    
    # Output results to show things off:
    rLogger.info("Start time: " + str(times[0]))
    rLogger.info("End time: " + str(times[-1]))
    rLogger.info("Radar | Fov (f or b), beam_start, gate_start, beam_end, gate_end")
    for r in relevant_radars:
        rLogger.info(r + ': ' + str(relevant_radars[r]))
    return conjs

def tag_conjunctions(fname):
    """
    Tags the ephemeris points in an RRI hdf5 file with SuperDARN 
    radars deemed to be nearby.

    *** HDF5 Tag Structure ***
    Added as a top-level group to the HDF5 file is a 
    'SuperDARN Conjunctions' group.

    Within this group will be groups for each radar whose FOV is 
    considered to be intersected by the CASSIOPE Ephemeris. The radars 
    will be listed by their SuperDARN Radar codes.

    Within these groups will be:
        - "Beam Extrema", representing bm_start and bm_end which are 
            the edge beams of the SuperDARN FOV which are intersected 
            in the dataset
        - "Gate Extrama", representing gt_start and gt_end which are 
            the extremes of the range gates of the SuperDARN FOV which 
            are intersected in the dataset.
        - "Intersection Quality", a 2 (good), 1 (moderate), or 0 (weak)
            which indicates the proximity and degree of intersection of
            the ephemeris through the FOV. (ideal: cross multiple beams
            and be within 40 range gates)
        - "Ephemeris Index Extrema", representing idx_start and idx_end
            which are the first and last indices in the ephemeris data
            of the hdf5 file that are formally intersecting the FOV.
        - "Radar Name", showing the full name of this radar
        - "Front or Back", showing whether or not the ephemeris track
            intersects the front or back of the radar FOV
    """
    import h5py
    # First, acquire the conjunctions from this RRI data file (get_conjs)
    conjs = get_conjunctions(fname)

    # Get their heuristic results
    hs = [ h.conj_heuristic() for h in conjs ] 

    # Create an 'SuperDARN Conjunctions' group in h5py and add subgroups
    # for each conjunction
    f = h5py.File(fname)
    f.create_group('SuperDARN Conjunctions')
    g = f['SuperDARN Conjunctions']

    # In each subgroup (titled by radar code), add 1-value datasets for
    # SuperDARN radar name, ephemeris start and end indices, bm_intersect_start
    # and bm_intersect_end, gt_intersect_start and gt_intersect_end
    for i, c in enumerate(conjs):
        g.create_group(str(c.code))
        gi = g[str(c.code)]
        gi.create_dataset('Radar Name', data=[str(c.name)])
        gi.create_dataset('Intersection Quality', data=[hs[i]])
        gi.create_dataset('Ephemeris Index Extrema', data=[c.idx_start, c.idx_end])
        gi.create_dataset('Beam Extrema', data=[c.bm_start, c.bm_end])
        gi.create_dataset('Gate Extrema', data=[c.gt_start, c.gt_end])
        gi.create_dataset('Front or Back', data=[c.forb])
    pass

def eliminate_conjunctions(fname):
    """
    Takes a filename for an RRI datafile which contains a 'SuperDARN 
    Conjunctions' data group, and deletes the datagroup.
    """
    import h5py
    f = h5py.File(fname)
    del f['SuperDARN Conjunctions']

def read_conjunctions(fname):
    """
    Takes an hdf5 file with SuperDARN Conjunctions tagged and extracts
    the conjunctions objects from it.
    """
    import h5py
    f = h5py.File(fname)
    try:
        g = f['SuperDARN Conjunctions']
    except KeyError:
        rLogger.error("No SuperDARN Conjunction tags in RRI File")
        return
    conjs = []
    for k in g.keys():
        ci = g[k]
        code = k
        name = ci['Radar Name'].value[0]
        bm_st, bm_end = ci['Beam Extrema'].value
        gt_st, gt_end = ci['Gate Extrema'].value
        idx_st, idx_end = ci['Ephemeris Index Extrema'].value
        qual = ci['Intersection Quality'].value[0]
        forb = ci['Front or Back'].value[0]
        conj = RRISuperdarnConjunction(code, name, forb, bm_st, gt_st,
                                       bm_end, gt_end, idx_st, idx_end)
        conjs.append(conj)
    return conjs

def fetch_radar_logs(fname):
    """
    Takes an RRI data file which has SuperDARN Conjunction tags added,
    fetches SuperDARN 'errlog' files for any Canadian radars that are 
    intersected, writing an output file that summarizes the conjunctions
    and the relevant lines of the errlog records.
    """
    nw = pydarn.radar.network()
    conjs = read_conjunctions(fname)
    uofs_rads = []
    for c in conjs:
        if c.name in ['Saskatoon', 'Prince George', 'Clyde River', 'Inuvik', 'Rankin Inlet']:
            uofs_rads.append(c)
  
    # Load data 
    data_path, data_fname = data_utils.initialize_data() 
    lons, lats, alts, ephtimes = data_utils.get_rri_ephemeris(data_fname)
    times = data_utils.ephems_to_datetime(ephtimes)
   
    # Start creating the format strings
    start = data_utils.ephem_to_datetime(ephtimes[0])
    end = data_utils.ephem_to_datetime(ephtimes[-1])
    start_str, end_str = start_end_strings(start, end)

    for u in uofs_rads:
        rcode = u.code
        # import rritk.plotting.plots as plots
        #plots.plot_fov_sat(u.name, lons, lats, date=start, suppress_show=True)

        f = data_utils.open_errlog(data_path, rcode, start) 

        #TODO: Consider actually copying this Errlog file
        start_line, end_line = start_end_lines(f, start_str, end_str)
        rLogger.info("Start line: {0}\nEnd line: {1}".format(start_line, end_line))

        # Having determined the relevant line of text in the errlog file, grab all the
        # relevant errlog data spanning the ephemeris file
        f.seekline(start_line) 
        rel_lines = f.readline()
        while f.line <= end_line:
            rel_lines = rel_lines + f.readline()

        time_tag = str(start.year) + two_pad(start.month) + \
            two_pad(start.day) +  "_" + two_pad(start.hour) + \
            two_pad(start.minute) 
        out_fname = OUTPUT_DIR + "/" + time_tag + "_" + rcode + ".dat"

        with open(out_fname, "w+") as f:
            f.write("OUTPUT FILE FOR RRI CONJUNCTION SCRIPT\n")
            f.write("======================================\n")
            f.write("Start Time: " + str(start) + "\nEnd Time: " + str(end))
            f.write("\nRRI File: "+ fname + "\n")
        
            for c in conjs:
                f.write("\n" + str(c))
            f.write("\n\nRelevant SuperDARN data from +/- 3 minutes of the RRI data:\n")
            f.write(rel_lines)
    
    return uofs_rads    

def start_end_strings(start, end):
    """
    Takes datetime objects for the start and end of an RRI dataset, returns
    the search strings to look for in SuperDARN .errlog files.
    """
    # ERRLOG files contain entries describing the beam number and frequency of 
    # the transmission, occurring roughly every three seconds, starting on each 
    # minute (at the ":00" mark).
    #
    # FURTHERMORE: the fact that the timestamps on SuperDARN data are untrustworthy
    # means we need to look at SuperDARN information up to a few minutes before and
    # after the RRI times of interest.
    #
    # Thus, this script will return all ERRLOG data starting *3* minutes before
    # the RRI's first timestamp, and ending *3* minutes after the RRI's last
    # timestamp, so as to give us some margin of error. The frequency that were
    # transmitted upon will hopefully inform the user as to where the two time
    # logs actually sync up.
    
    # Also, need to pad the datetime objects so that their format matches the 
    # filename and timestamp formats for the ERRLOG files.
    st_month = "0" + str(start.month) if str(start.month).__len__() == 1 else str(start.month)
    st_day = "0" + str(start.day) if str(start.day).__len__() == 1 else str(start.day)
    st_hour = "0" + str(start.hour) if str(start.hour).__len__() == 1 else str(start.hour)
    end_hour = "0" + str(end.hour) if str(end.hour).__len__() == 1 else str(end.hour)
    # Also, due to untrustworthy timestamps, we look at times in the ERRLOG file
    # up to 3 minutes before and after the RRI times.
    st_min = start.minute-3
    end_min = end.minute+3
    st_min = "0" + str(st_min) if str(st_min).__len__() == 1 else str(st_min)
    end_min = "0" + str(end_min) if str(end_min).__len__() == 1 else str(end_min)
    
    start_string = " " + str(st_hour) + ":" +  str(st_min) + ":" # "00" # Omit seconds for now in case
    end_string = " " + str(end_hour) + ":" + str(end_min) + ":" #"00" # they don't follow expected pattern
    rLogger.info("Start and End strings to search for: {0}, {1}".format(start_string, end_string))
    return start_string, end_string

def start_end_lines(f, start_string, end_string):
    """
    Searches the errlog file pointed to by file handle f for the 
    starting and ending intervals
    """
    # With the file open, search for the desired time interval's beginning.
    # Search for the position of the first line of interest
    found = False
    while not found:
        ln = f.readline()
        #rLogger.debug("Line just read: {0}".format(ln))
        if ln.find(start_string) != -1:
            found = True
            start_line = f.line
            rLogger.info(str(start_line) + ": " + ln)
    # Now search for the position of the final line of interest
    found = False
    while not found:
        ln = f.readline()
        if ln.find(end_string) != -1:
            found = True
            end_line = f.line
            rLogger.info(str(end_line) + ": " + ln)
        elif f.line > (start_line + 400): #1000000 bytes
            end_line = start_line + 400 
            found = True
            rLogger.info(str(end_line) + ": " + ln)

    return start_line, end_line

# ----------------------------------------------------------------------------- 

def initialize_logger(quiet_mode=False):
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

if __name__ == "__main__":
    """
    If script is run directly with a RRI .h5 filename as an argument, then
    the desired execution is to load the ephemeris points and calculate 
    intersections and spit them out.
    """
    import sys
    initialize_logger()
    rLogger.info(sys.argv)
    if len(sys.argv) > 1 and type(sys.argv[1]) == str:
    # An additional argument has been provided
        try:
            rLogger.info("Attempting to read from .h5 file provided as argument...")
            lons, lats, alts, ephtimes = data_utils.get_rri_ephemeris(sys.argv[1])
        except:
            rLogger.info("Catch-all except statement. Something didn't work with loading the provided arg file")
            exit()
    try:
        rLogger.info("Attempting to get conjunctions then!")
        conjs = get_conjunctions(sys.argv[1])
        rLogger.info("Conjunction output looks like: {0}".format(conjs))
    except Exception as e:
        rLogger.error("Getting conjunctions didn't work so well, exception thrown: {0}".format(e))
        pass
