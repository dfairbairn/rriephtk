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

import data_utils
import script_utils
import range_cells

logging.basicConfig(filename='conjunctions.log',level=logging.ERROR)

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
    times = script_utils.ephems_to_datetime(ephtimes)

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
                start = non_nan_b[0] 
                end = non_nan_b[ non_nan_b.__len__() - 1 ]
                conjs.append(RRISuperdarnConjunction(
                    rad.code[0], rad.name, 'back', bm_b[start], gt_b[start], 
                    bm_b[end], gt_b[end], idx_start=start, idx_end=end))
                relevant_radars[rad.name] = (rad.code[0],"back", bm_b[start],gt_b[start],bm_b[end],gt_b[end])

            if (non_nan_f.__len__() > 0):
                start = non_nan_f[0]
                end = non_nan_f[ non_nan_f.__len__() - 1 ]
                conjs.append(RRISuperdarnConjunction(
                    rad.code[0], rad.name, 'front',  bm_f[start], gt_f[start], 
                    bm_f[end], gt_f[end], idx_start=start, idx_end=end))
                relevant_radars[rad.name] = (rad.code[0],"front",bm_f[start],gt_f[start],bm_f[end],gt_f[end])
        
        #Show progress on screen
        sys.stdout.flush()
        script_utils.update_progress((i+1)/float(len(activerads)))

    end_t = timeit.default_timer()
    print("\nTime req'd to compute detailed intersections by brute force: ",
        str(end_t - start_t), " seconds.")

    # In general there will be a post-processing step in which only the 
    # conjunction results will be included in the output
    
    # Output results to show things off:
    print("Start time: " + str(times[0]))
    print("End time: " + str(times[-1]))
    print("Radar | Fov (f or b), beam_start, gate_start, beam_end, gate_end")
    for r in relevant_radars:
        print r + ': ' + str(relevant_radars[r])
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
    return None

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
        logging.error("No SuperDARN Conjunction tags in RRI File")
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
    intersected.
    """

    conjs = read_conjunctions(fname)
    uofs_rads = []
    for c in conjs:
        if c.name in ['Saskatoon', 'Prince George', 'Clyde River', 'Inuvik', 'Rankin Inlet']:
            uofs_rads.append(c)
   
    data_path, data_fname = data_utils.initialize_data() 





    return uofs_rads    

if __name__ == "__main__":
    """
    If script is run directly with a RRI .h5 filename as an argument, then
    the desired execution is to load the ephemeris points and calculate 
    intersections and spit them out.
    """
    import sys
    print(sys.argv)
    if len(sys.argv) > 1:
        # An additional argument has been provided 
        if type(sys.argv[1])==str:
            try:
                print("Privyet")
                lons, lats, alts, ephtimes = data_utils.get_rri_ephemeris(sys.argv[1])
            except:
                print("Catch-all except statement. Something didn't work with loading the provided arg file")
                exit()

    conjs = get_conjunctions(sys.argv[1])
 
