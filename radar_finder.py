"""
file: 'radarFinder.py'

description:
  Contains functions for determining which SuperDARN radar fields-of-view 
  contain a given input geographic point. 

date: June 3rd, 2016

author: David Fairbairn


""" """
--------------------------------- PART IV --------------------------------------
                                RADAR FINDER         
--------------------------------------------------------------------------------
"""
import numpy as np
import davitpy
from datetime import *

import logging

def radar_finder(lats, lons):
    """
    Given latitude and longitude arrays, returns the SuperDARN radars whose front
    or back lobes contain the coordinates.
    
    -- PARAMS --
    lats: list or array of input latitude points
    
    lons: list or array of input longitude points
    
    
    -- RETURNS -- 
    results: a list of the same length as the lats/lons lists, which stores lists
        containing the names of the relevant SuperDARN radars for that coordinate.
    
    
    written by David Fairbairn, 2016-05
    """
    assert(np.shape(lats) == np.shape(lons)),"Lat/Long arrays must be same length."
    
    lons,lats = (np.array(lons)+360.)%360.0,np.array(lats)
    num = np.size(lons)
    
    nw = davitpy.pydarn.radar.network()
    results = dict()
    for i in range(np.size(lons)):
        for rad in nw.radars:
            if rad.status == 1:
                tst_containment = in_fov(lats[i],lons[i],rad_id=rad.id)




def in_fov(lat, lon, rad_code=None,rad_id=None):
    """
    Checks if the lat/lon coordinate pair is within the given radar's front fov
    or back lobe.

    Similar to Angeline Burrell's pos_to_rg() function, it requires either for
    rad_code or rad_id to be supplied.
    
    -- PARAMS --
    lat: (double)
        For now, a *single* latitude value in degrees
    lon: (double)
        For now, a *single* longitude value in degrees
    rad_code: (String or NoneType)
        The 3-letter string denoting a radar. Either this or rad_id is required.
    rad_id: (int or NoneType)
        The integer id for a radar. Either this or rad_code is required.

    -- RETURNS -- 
    in_fov: (boolean)
        Boolean marker of whether or not the coordinate lies in the radar's range 
    

    """

#    assert isinstance(rad_code, str) or isinstance(rad_id, int), \
#        logging.error('must specify either rad_code or rad_id')

    # if both are provided, default to using the code
#    if isinstance(rad_code, str) and isinstance(rad_id, int):
#        rad_id = None

    # ensure longitude is in the same format as used in DaVitPy
    lon = (lon+360.)%360

    # Grab the site details (gets antenna lat/lon)
    hw = davitpy.pydarn.radar.site(code=rad_code, radId=rad_id, dt=datetime.now()) 
    sides = ['front','back']

    containment = []
    for f in sides:
        # Grab FOV via the radar structure
        fov = davitpy.pydarn.radar.radFov.fov(site=hw,altitude=300.,model='IS',
                                        coords='geo',ngates=75,fov_dir=f)
    
        # Convert query point
        leftFovLat = fov.latFull[0][75]
        leftFovLon = fov.lonFull[0][75]
    
        rightFovLat = fov.latFull[16][75]
        rightFovLon = fov.lonFull[16][75]
    
        centerlat = hw.geolat
        centerlon = (hw.geolon+360.)%360.
        centerHemisphere = 1 if (centerlat >= 0) else -1
    
        queryHemisphere = 1 if (lat >= 0) else -1   
        queryX, queryY = get_xy(lat, lon)

        # no radar is considered to be capable of reaching the other hemisphere 
        if centerHemisphere != queryHemisphere:
            return False 
    
        centerX, centerY = get_xy(centerlat, centerlon)
        leftFovX, leftFovY = get_xy(leftFovLat, leftFovLon)
        rightFovX, rightFovY = get_xy(rightFovLat, rightFovLon)

        leftDist = np.linalg.norm(np.array((leftFovX,leftFovY)) - np.array((centerX,centerY)))   
        rightDist = np.linalg.norm(np.array((rightFovX,rightFovY)) - np.array((centerX,centerY)))
        pointDist = np.linalg.norm(np.array((queryX,queryY)) - np.array((centerX,centerY)))

        if pointDist > max(leftDist,rightDist):
            return False # query point can't be farther than FOV's range away from center

        # Now that we have everything in a flat plane, begin testing conditions.
    
        # to be in radar front fov, the query point should be RIGHT of the left line
        # and to the LEFT. Back FOV: LEFT of left line, RIGHT of right. 
       
        # Note: if the left boundary line is positive, this would mean wanting to
        # be 'under' it, whereas if  it were negative we would want to be above.
        # Then we would have to consider the right boundary as well for each case,
        # and then these conditions are reversed for the back lobe. So it's instead
        # more practical to simply calculate the relative angles between the query
        # point and the left and right boundaries.
     
        #Put the tests for a long/lat being in a FOV into a generalized function
    
        containment.append((f,point_in_sector(queryX,queryY,centerX,centerY,leftFovX,leftFovY,rightFovX,rightFovY)))

    return containment



def get_xy(lat, lon): 
    """
    Takes a latitude and longitude coordinate and converts it to x & y 
    coordinates (in a sphere's equatorial plane).

    -- PARAMS --
    lat: latitude on a sphere
    lon: longitude on a sphere
    """
    lon = (lon+360.)%360
    x = (90 - np.abs(lat)) * np.sin(np.deg2rad(lon))
    y = -(90 - np.abs(lat)) * np.cos(np.deg2rad(lon))
    return (x,y)

def point_in_sector(pointX,pointY,centerX,centerY,leftX,leftY,rightX,rightY,verbosity="off"):
    """
    General function for checking if a point is in a certain sector of a circle
    defined around that point in the plane.

    -- PARAMS --
    pointX, pointY: a point that we're testing
    centerX, centerY: the point at the center of the circle
    leftX, leftY: a point along the left edge of the sector.
    rightX, rightY: a point along the right edge of the sector.

    -- RETURNS -- 


    """
    leftDist = np.linalg.norm(np.array((leftX,leftY)) - np.array((centerX,centerY)))   
    rightDist = np.linalg.norm(np.array((rightX,rightY)) - np.array((centerX,centerY)))
    pointDist = np.linalg.norm(np.array((pointX,pointY)) - np.array((centerX,centerY)))
 
#    leftDist = ((leftX - centerX)**2.0 + (leftY - centerY)**2.0)**0.5
#    rightDist = ((rightX - centerX)**2.0 + (rightY - centerY)**2.0)**0.5
#    pointDist = ((pointX - centerX)**2.0 + (pointY - centerY)**2.0)**0.5

    if pointDist > max(leftDist,rightDist):
        return False # query point can't be farther than FOV's range away from center

    # Use Law of Cosines to determine angles: angle C = arccos( (a^2 + b^2 - c^2)/(2ab) )
    distFromLeft = np.linalg.norm(np.array((leftX,leftY)) - np.array((pointX,pointY)))
    angleFromLeft = np.arccos((leftDist**2 + pointDist**2 - distFromLeft**2)/(2*leftDist*pointDist))
  
    distFromRight = np.linalg.norm(np.array((rightX,rightY)) - np.array((pointX,pointY)))
    angleFromRight = np.arccos((rightDist**2 + pointDist**2 - distFromRight**2)/(2*rightDist*pointDist)) 

    distLeftRight = np.linalg.norm(np.array((leftX,leftY)) - np.array((rightX,rightY)))
    angleLeftRight = np.arccos((rightDist**2 + leftDist**2 - distLeftRight**2)/(2*leftDist*rightDist))

    #if (verbosity != "off"):
    #    verdict = abs(angleLeftRight - angleFromLeft - angleFromRight) < 0.001
    #    print "Angle from left edge of sector to line: " + str(angleFromLeft)
    #    print "Angle from right edge of sector to line: " + str(angleFromRight)
    #    print "Angle between left and right edges of sector: " + str(angleLeftRight)
    #    print "Verdict: " + str(verdict)

    return (abs(angleLeftRight - angleFromLeft - angleFromRight) < 0.001)








"""
TESTING
-------





"""
if __name__ == "__main__":

    # TESTING get_xy()
    



    # TESTING line_in_sector()
    n = 10
    N = 2*n + 1
    list_compr = [(x,y) for x in (np.array(range(N)) - n) for y in (np.array(range(N)) - n)]
    successes = []
    for p in list_compr:
        if point_in_sector(p[0],p[1],0,0,4,4,4,3):
            successes.append(p) 


    # TESTING in_fov()

    #lst = [(x,y) for x in (np.array(range(5))/10.00 + 52.13) for y in (np.array(range(5))/10.00 + 256.67)]
    #for p in lst:
    #    results = in_fov(p[0],p[1],rad_code='sas')
    #print "Testing 5000 samples took: " + str(end - start) + " seconds."
    
    import timeit
    start = timeit.default_timer()
    in_fov(54.0,257.0,rad_code='sas')
    end = timeit.default_timer()

    # TESTING radar_finder()
    print "Testing the full radar finder with a subset of RRI Ephemeris latitude/longitude pairs..."
    import h5py
    dat_fname = "/home/david/pyth_stuff/script/data/RRI_20160401_072714_073111_lv1_v1.h5" 
    f = h5py.File(dat_fname)
    
    geog_longs = f['CASSIOPE Ephemeris']['Geographic Longitude (deg)'].value
    geog_lats = f['CASSIOPE Ephemeris']['Geographic Latitude (deg)'].value

    lons = geog_longs[0:10]
    lats = geog_lats[0:10]

    # 210 seconds when distance calculations and tests are in point_in_sector()
    # 171 seconds when distance test is before calling point_in_sector (but req'ing more calcs)
    # TODO:
    # ??? seconds by getting rid of calls to get_xy() ???
    # ??? seconds by getting rid of point_in_sector() and just having everything in one procedure??


    start = timeit.default_timer()
    results = radar_finder(lats,lons)
    end = timeit.default_timer()

    print "Computing relevant radars for 11 Ephemeris points took : " + str(end-start) + " seconds."
