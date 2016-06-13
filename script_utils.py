"""
File: 'script-utils.py'
Description: Contains utilities used in the RRI Conjunction Finder script.
    These include functions to more conveniently grab FOV data from Davitpy,
    plotting functions that make use of Matplotlib, etc.
Author: David Fairbairn
Date: June 13th, 2016


"""

def plot_fov_sat(fov, ephem_lons, ephem_lats):
    """
    



    """
    import matplotlib
    fovlons = ((fov.lonFull+360.)%360.).ravel()
    fovlats = (fov.latFull).ravel()
    fovcol = 3*np.ones(np.size(fovlons))

    ephemcol = 5*np.ones(np.size(ephem_lats))
    ephemlons = (ephem_lons+360.)%360.
    ephemlats = ephem_lats

    colors = np.concatenate((fovcol,ephemcol))
    lons = np.concatenate((fovlons,ephemlons))
    lats = np.concatenate((fovlats,ephemlats))

    matplotlib.pyplot.scatter(lons,lats,c=colors,edgecolors='face')
    matplotlib.pyplot.show()

def plot_fovs_sat(fovs, ephem_longs, ephem_lats):



    import matplotlib
    if (np.shape(fovs)[0] > 10):
        print "Can't do more than 10 FOVs"
        return

    lons = []
    lats = []
    colors = []
    color_offset = 0
    for fov in fovs:
        fovlons = ((fov.lonFull+360.)%360.).ravel()
        fovlats = (fov.latFull).ravel()
        fovcol = (1+color_offset)*np.ones(np.size(fovlons))
        color_offset += 1
        
        lons = np.concatenate((lons,fovlons))
        lats = np.concatenate((lats,fovlats))
        colors = np.concatenate((colors,fovcol))
    
    ephemlons = (ephem_longs+360.)%360.
    ephemlats = ephem_lats 
    ephemcol = (color_offset+1)*np.ones(np.size(ephemlons))
    lons = np.concatenate((lons,ephemlons)) 
    lats = np.concatenate((lats,ephemlats))
    colors = np.concatenate((colors,ephemcol))
    
    matplotlib.pyplot.scatter(lons,lats,c=colors,edgecolors='face')
    matplotlib.pyplot.show()


def get_fov_by_name(name):
    from davitpy import pydarn
    nw = pydarn.radar.network()
    rad_id = (nw.getRadarByName(name)).id
    site = pydarn.radar.site(radId = rad_id)
    fov = pydarn.radar.radFov.fov(site=site,altitude=300.0,model='IS',coords='geo',ngates=75)
    return fov


