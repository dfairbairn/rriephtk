#!/usr/bin/env python
# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------------
# rgCoords.py, Angeline G. Burrell (AGB), UoL
#
# Comments: Scripts to convert data to and from range gate/beam coordinates
#-----------------------------------------------------------------------------
"""convert data to and from range gate/beam coordinates

Functions
------------------------------------------------------------------------------
pos_to_rg   Convert from geographic or magnetic coordinates to range gate/beam
------------------------------------------------------------------------------

Author: Angeline G. Burrell (AGB)
Date: February 25, 2016
Inst: University of Leicester (UoL)

"""
import logging
import numpy as np
import datetime as dt
import time

def pos_to_rg(lat, lon, alt, coords, dtime, frang, rsep, bmsep=None,
              maxgate=None, maxbeam=None, rad_code=None, rad_id=None,
              ret_fov="both"):
    """Convert from geographic or magnetic coordinates to range gate and beam

    Parameters
    --------------------
    lat : (np.array)
        Array of latitudes in degrees
    lon : (np.array)
        Array of longitudes in degrees or hours
    coords : (str)
        String denoting the coordinates.  Choose from "geo", "mag", and "mlt"
    alt : (float)
        Altitude of data in km above the surface of the earth
    dtime : (dt.datetime)
        Time at which to calculate positions
    frang : (float)
        Distance to first range gate (km)
    rsep : (float)
        Seperation between range gates (km)
    bmsep : (float)
        Seperation between beams.  If None, use value from radar hardware file.
        (default=None,)
    maxgate : (int or NoneType)
        Maximum number of gates.  If None, use value from radar hardware file.
        (default=None)
    maxbeam : (int or NoneType)
        Maximum number of beams.  If None, use value from radar hardware file.
        (default=None)
    rad_code : (string or NoneType)
        radar id as a 3 character string.  Must enter either rad_code or
        rad_id.  If both are input, rad_code will be used. (default=None)
    rad_id : (int or NoneType)
        radar id an interger.  Must enter either rad_code or
        rad_id.  If both are input, rad_code will be used. (default=None)
    ret_fov : (str)
        string denoting which radar field-of-view to return ('front', 'back',
        or 'both').  (default='both')

    Returns
    -----------
    range_gate : (np.array)
        Array of range gates
    bmnum : (np.array)
        Array of beam numbers
    fovflg : (np.array)
        Array of field-of-view flags, where 1 indicates front and -1 indicates
        back.
    """
    import davitpy.pydarn.radar as pyrad
    import davitpy.utils.geoPack as geo
    tt = time.clock() # TIME

    #------------------
    # Test input
    assert isinstance(lat, np.ndarray), \
        logging.error('lat must be a numpy array')
    assert isinstance(lon, np.ndarray), \
        logging.error('lon must be a numpy array')

    if lat.shape != lon.shape:
        logging.error('latitude and longitude arrays must be the same size')
        return None, None, None

    assert isinstance(rad_code, str) or isinstance(rad_id, int), \
        logging.error('must specify either rad_code or rad_id')

    if isinstance(rad_code, str) and isinstance(rad_id, int):
        rad_id = None

    print "TIME: test input", time.clock() - tt
    tt = time.clock()

    #---------------------------------------------------
    # Initialize the output
    range_gate = np.ones(shape=lat.shape, dtype=int) * -1
    bmnum = np.ones(shape=lat.shape, dtype=int) * -1
    fovflg = np.zeros(shape=lat.shape, dtype=int)

    print "TIME: initialize output", time.clock() - tt
    tt = time.clock()
    #----------------------------------
    # Initialize radar information
    hard = pyrad.site(code=rad_code, radId=rad_id, dt=dtime)
    (d1, d2, radius) = geo.geodToGeoc(hard.geolat, hard.geolon, inverse=True)

    if maxbeam is None:
        maxbeam = hard.maxbeam

    if maxgate is None:
        maxgate = hard.maxgate

    if bmsep is None:
        bmsep = hard.bmsep

    ivertices = [[1, 2], [0, 3], [0, 3], [1, 2]]

    print "TIME: initialize radar", time.clock() - tt
    tt = time.clock()

    #----------------------------------
    # Get the radar field-of-views
    fov_val = {"front":1, "back":-1}
    get_fov = [f for f in fov_val.keys() if f == ret_fov or ret_fov == "both"]
    limits = np.ones(shape=(maxbeam, 1+maxgate*2, 4, 2), dtype=float) * np.nan
    center = np.ones(shape=(maxbeam, 1+maxgate*2, 2), dtype=float) * np.nan
    cdist = np.ones(shape=(maxbeam, 1+maxgate*2), dtype=float) * np.nan

    for f in get_fov:
        rad_fov = pyrad.radFov.fov(site=hard, frang=frang, rsep=rsep,
                                   nbeams=maxbeam, ngates=maxgate, bmsep=bmsep,
                                   altitude=alt, model="IS", coords=coords,
                                   date_time=dtime, fov_dir=f)

        for b in range(maxbeam):
            for r in range(maxgate):
                # Index for range gate isn't the range gate because indexes
                # can't be negative and need to account for field-of-view
                rg = maxgate - r - 1 if f == "back" else r + maxgate

    
                # Save the location of the central point
                center[b][rg][0] = rad_fov.latCenter[b][r]
                center[b][rg][1] = rad_fov.lonCenter[b][r]

                # Save the vertex coordinates 
                # beam min, range gate min
                limits[b][rg][0][0] = rad_fov.latFull[b][r]
                limits[b][rg][0][1] = rad_fov.lonFull[b][r]
                d1 = geo.greatCircleDist(center[b][rg][0], center[b][rg][1],
                                         limits[b][rg][0][0],
                                         limits[b][rg][0][1]) * radius
                # beam min, range gate max
                limits[b][rg][1][0] = rad_fov.latFull[b][r+1]
                limits[b][rg][1][1] = rad_fov.lonFull[b][r+1]
                d2 = geo.greatCircleDist(center[b][rg][0], center[b][rg][1],
                                         limits[b][rg][1][0],
                                         limits[b][rg][1][1]) * radius
                # beam max, range gate min
                limits[b][rg][2][0] = rad_fov.latFull[b+1][r]
                limits[b][rg][2][1] = rad_fov.lonFull[b+1][r]
                d3 = geo.greatCircleDist(center[b][rg][0], center[b][rg][1],
                                         limits[b][rg][2][0],
                                         limits[b][rg][2][1]) * radius
                # beam max, range gate max
                limits[b][rg][3][0] = rad_fov.latFull[b+1][r+1]
                limits[b][rg][3][1] = rad_fov.lonFull[b+1][r+1]
                d4 = geo.greatCircleDist(center[b][rg][0], center[b][rg][1],
                                         limits[b][rg][3][0],
                                         limits[b][rg][3][1]) * radius

                # Save the largest distance between a vertex and the center
                cdist[b][rg] = max([d1, d2, d3, d4])

    print "TIME: set fov info", time.clock() - tt
    tt = time.clock()
    #-----------------------------------------------------------------------
    # Cycle through the locations, determining which points fall within the
    # radar fields-of-view and where they lie
    for i,l in enumerate(lat):
        # Cycle through all possible bins, testing to see if 
        for b in range(limits.shape[0]):
            for r in range(limits.shape[1]):
                if not np.any(np.isnan(center[b][r])):
                    # Test to see if the distance between this point and the
                    # center is at least as small as than the maximum distance
                    # between a vertex and the center
                    d = geo.greatCircleDist(center[b][r][0], center[b][r][1],
                                            l, lon[i]) * radius

                    if d <= cdist[b][r]:
                        # Other possible indices with acceptable distances
                        ib = [b, b-1, b-1, b-1, b, b, b+1, b+1, b+1]
                        ir = [r, r-1, r, r+1, r-1, r+1, r-1, r, r+1]
                        dd = [np.nan for rg in ir]
                        ii = [True for rg in ir]
                        dd[0] = d

                        for j in np.arange(1, 9):
                            rg = ir[j]
                            bm = ib[j]
                            if(rg >= 0 and rg < limits.shape[1] and
                               bm >= 0 and bm < limits.shape[0] and
                               not np.any(np.isnan(center[bm][rg]))):
                                d1 = geo.greatCircleDist(center[bm][rg][0],
                                                         center[bm][rg][1],
                                                         l, lon[i]) * radius
                                if d1 <= cdist[bm][rg]:
                                    dd[j] = d1
                            else:
                                ii[j] = False

                        d1 = np.nanmin(dd)
                        d2 = dd.index(d1)
                        rg = ir[d2]
                        bm = ib[d2]

                        # Now that we've found the shortest distance between
                        # a range gate center and the specified point, and we
                        # know that this distance is small enough to fit in
                        # this range gate-beam box, see whether it actually
                        # does.
                        d3 = True if np.all(ii) else False

                        if not d3:
                            # At least one edge of this box is a field-of-view
                            # boundary.  Test to see of the specified point
                            # falls outside of this boundary.
                            #
                            # The first step is to identify the box vertex
                            # closest to the specified point and express it
                            # in polar coordinates
                            dd = [geo.greatCircleDist(ll[0], ll[1], l, lon[i])
                                  * radius for ll in limits[bm][rg]]
                            iv = dd.index(min(dd))
                            c2_v = [np.cos(np.radians(limits[bm][rg][iv][0])) *
                                    np.sin(np.radians(limits[bm][rg][iv][1])),
                                    np.cos(np.radians(limits[bm][rg][iv][0])) *
                                    np.cos(np.radians(limits[bm][rg][iv][1])),
                                    np.sin(np.radians(limits[bm][rg][iv][0]))]

                            # Locate the intercepts between the great circles
                            # that pass through the central and specified points
                            # and the walls of the box that pass through the
                            # closest vertex
                            d4 = np.nan # set the minimum distance to np.nan

                            # Define the great circle between the center and
                            # the specified point, starting with expressing
                            # the central and specified point in polar coords.
                            c1_p1 = [np.cos(np.radians(center[bm][rg][0])) *
                                     np.sin(np.radians(center[bm][rg][1])),
                                     np.cos(np.radians(center[bm][rg][0])) *
                                     np.cos(np.radians(center[bm][rg][1])),
                                     np.sin(np.radians(center[bm][rg][0]))]
                            c1_p2 = [np.cos(np.radians(l)) *
                                     np.sin(np.radians(lon[i])),
                                     np.cos(np.radians(l)) *
                                     np.cos(np.radians(lon[i])),
                                     np.sin(np.radians(l))]

                            # Get normal to planes containing this great circle
                            # by taking the cross product of the vector to each
                            # point from the origin.
                            n1 = np.cross(c1_p1, c1_p2)

                            for iw in ivertices[iv]:
                                wlat = limits[bm][rg][iw][0]
                                wlon = limits[bm][rg][iw][1]

                                # Define the great circle between the selected
                                # and adjoining vertices, starting by expressing
                                # the adjoining vertex in polar coordinates
                                c2_w = \
                                    [np.cos(np.radians(limits[bm][rg][iw][0])) *
                                     np.sin(np.radians(limits[bm][rg][iw][1])),
                                     np.cos(np.radians(limits[bm][rg][iw][0])) *
                                     np.cos(np.radians(limits[bm][rg][iw][1])),
                                     np.sin(np.radians(limits[bm][rg][iw][0]))]

                                # Get the normal plane for this great circle
                                n2 = np.cross(c2_v, c2_w)

                                # Find the line of intersection between two
                                # planes, which is the normal to the poles of
                                # each plane.
                                line = np.cross(n1, n2)

                                # Find intersection point in polar coordinates,
                                # the other intersection is at the negative
                                # location of this point
                                pint = line / np.sqrt(np.dot(line, line))

                                # Convert from polar coordinates to latitude
                                # and longitude
                                ilat1 = np.degrees(np.arcsin(pint[2]))
                                ilon1 = np.degrees(np.arctan2(pint[1], pint[0]))

                                pint *= -1.0
                                ilat2 = np.degrees(np.arcsin(pint[2]))
                                ilon2 = np.degrees(np.arctan2(pint[1], pint[0]))

                                di1 = geo.greatCircleDist(center[bm][rg][0],
                                                          center[bm][rg][1],
                                                          ilat1, ilon1) * radius
                                di2 = geo.greatCircleDist(center[bm][rg][0],
                                                          center[bm][rg][1],
                                                          ilat2, ilon2) * radius

                                # Save only the minimum distance
                                if di2 < di1:
                                    di1 = di2
                                if np.isnan(d4) or di1 < d4:
                                    d4 = di1

                            # Compare the distance between the center and the
                            # closest intercept to the distance between the
                            # center and the specified point.  If the second
                            # distance is smaller, the specified point lies
                            # inside the box.
                            if d4 >= d1:
                                d3 = True

                        # Assign beam, range gate, and field-of-view flags
                        # if this point was found to lie inside the box
                        if d3:
                            bmnum[i] = bm
                            if rg <= maxgate:
                                fovflg[i] = -1
                                range_gate[i] = maxgate - rg - 1
                            else:
                                fovflg[i] = 1
                                range_gate[i] = rg - maxgate

    print "TIME: done", time.clock() - tt
    return bmnum, range_gate, fovflg

def great_circle_intersections(c1_lat1, c1_lon1, c1_lat2, c1_lon2, c2_lat1,
                               c2_lon1, c2_lat2, c2_lon2):
    """Find the intersection points for two great circles

    Parameters
    -----------
    c1_lat1 : (float)
        Latitude in degrees of one point that defines the first great circle
    c1_lon1 : (float)
        Longitude in degrees of one point that defines the first great circle
    c1_lat2 : (float)
        Latitude in degrees of one point that defines the first great circle
    c1_lon2 : (float)
        Longitude in degrees of one point that defines the first great circle
    c2_lat1 : (float)
        Latitude in degrees of one point that defines the first great circle
    c2_lon1 : (float)
        Longitude in degrees of one point that defines the first great circle
    c2_lat2 : (float)
        Latitude in degrees of one point that defines the first great circle
    c2_lon2 : (float)
        Longitude in degrees of one point that defines the first great circle

    Returns
    ---------
    ilat1 : (float)
        Latitude of first intersection point in degrees
    ilon1 : (float)
        Longitude of first intersection point in degrees
    ilat2 : (float)
        Latitude of second intersection point in degrees
    ilon2 : (float)
        Longitude of second intersection point in degrees

    Notes
    ------
    Each great circle defined by two points on a sphere.  If both great circles
    are identical, np.nan will be returned for all intersection points
    """
    # Express each point in polar coordinates [x, y, z]
    c1_p1 = [np.cos(np.radians(c1_lat1)) * np.sin(np.radians(c1_lon1)),
             np.cos(np.radians(c1_lat1)) * np.cos(np.radians(c1_lon1)),
             np.sin(np.radians(c1_lat1))]

    c1_p2 = [np.cos(np.radians(c1_lat2)) * np.sin(np.radians(c1_lon2)),
             np.cos(np.radians(c1_lat2)) * np.cos(np.radians(c1_lon2)),
             np.sin(np.radians(c1_lat2))]

    c2_p1 = [np.cos(np.radians(c2_lat1)) * np.sin(np.radians(c2_lon1)),
             np.cos(np.radians(c2_lat1)) * np.cos(np.radians(c2_lon1)),
             np.sin(np.radians(c2_lat1))]

    c2_p2 = [np.cos(np.radians(c2_lat2)) * np.sin(np.radians(c2_lon2)),
             np.cos(np.radians(c2_lat2)) * np.cos(np.radians(c2_lon2)),
             np.sin(np.radians(c2_lat2))]

    # Get normal to planes containing great circles by taking the cross product
    # of the vector to each point from the origin.
    n1 = np.cross(c1_p1, c1_p2)
    n2 = np.cross(c2_p1, c2_p2)

    # Find the line of intersection between two planes, which is the normal to
    # the poles of each plane.
    line = np.cross(n1, n2)

    # Find intersection point in polar coordinates (other intersection is
    # at the negative location of this point
    pint = line / np.sqrt(np.dot(line, line))

    # Convert from polar coordinates to latitude and longitude
    ilat1 = np.degrees(np.arcsin(pint[2]))
    ilon1 = np.degrees(np.arctan2(pint[1], pint[0]))

    pint *= -1.0
    ilat2 = np.degrees(np.arcsin(pint[2]))
    ilon2 = np.degrees(np.arctan2(pint[1], pint[0]))

    return ilat1, ilon1, ilat2, ilon2
