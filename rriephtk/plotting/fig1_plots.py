from rriephtk.utils.data_utils import *
from rriephtk.analysis.analysis_tools import *
from davitpy.utils import plotUtils

from mpl_toolkits import basemap

MILLSTONE_TX_LON = -71.491 # acquired from google maps satellite imagery
MILLSTONE_TX_LAT = 42.619

class mapObjPlus(basemap.Basemap):
  """This class wraps arround :class:`mpl_toolkits.basemap.Basemap` (<http://tinyurl.com/d4rzmfo>)
  
  **Members**: 
    * **coords** (str): map coordinate system ('geo', 'mag', 'mlt').
    * all members of :class:`mpl_toolkits.basemap.Basemap` (<http://tinyurl.com/d4rzmfo>) 
  **Methods**:
    * all methods of :class:`mpl_toolkits.basemap.Basemap` (<http://tinyurl.com/d4rzmfo>)
  **Example**:
    ::

      # Create the map
      myMap = utils.mapObj(boundinglat=30, coords='mag')
      # Plot the geographic and geomagnetic North Poles
      # First convert from lat/lon to map projection coordinates...
      x, y = myMap(0., 90., coords='geo')
      # ...then plot
      myMap.scatter(x, y, zorder=2, color='r')
      # Convert to map projection...
      x, y = myMap(0., 90., coords='mag')
      # ...and plot
      myMap.scatter(x, y, zorder=2, color='g')

  .. note:: Once the map is created, all plotting calls will be assumed to already be in the map's declared coordinate system given by **coords**.

  """

  def __init__(self, datetime=None, coords='geo', 
    projection='stere', resolution='c', dateTime=None, 
    lat_0=None, lon_0=None, boundinglat=None, width=None, height=None, 
    fillContinents='.8', fillOceans='None', fillLakes=None, coastLineWidth=0., 
    grid=True, gridLabels=True, showCoords=True, **kwargs):
    """Create empty map 
    
    **Args**:    
      * **[width, height]**: width and height in m from the (lat_0, lon_0) center
      * **[lon_0]**: center meridian (default is -70E)    
      * **[lat_0]**: center latitude (default is -90E)
      * **[boundingLat]**: bounding latitude (default it +/-20)    
      * **[grid]**: show/hide parallels and meridians grid    
      * **[fill_continents]**: continent color. Default is 'grey'    
      * **[fill_water]**: water color. Default is 'None'    
      * **[coords]**: 'geo'
      * **[showCoords]**: display coordinate system name in upper right corner
      * **[dateTime]** (datetime.datetime): necessary for MLT plots if you want the continents to be plotted
      * **[kwargs]**: See <http://tinyurl.com/d4rzmfo> for more keywords
    **Returns**:
      * **map**: a Basemap object (<http://tinyurl.com/d4rzmfo>)
    **Example**:
      ::

        myMap = mapObj(lat_0=50, lon_0=-95, width=111e3*60, height=111e3*60)
        
    written by Sebastien, 2013-02
    """
    from davitpy.models import aacgm
    from pylab import text
    import math
    from copy import deepcopy

    self._coordsDict = {'mag': 'AACGM',
              'geo': 'Geographic',
              'mlt': 'MLT'}

    if coords is 'mlt':             
      print('MLT coordinates not implemented yet.')
      return

    # Add an extra member to the Basemap class
    if coords is not None and coords not in self._coordsDict:
      print('Invalid coordinate system given in coords ({}): setting "geo"'.format(coords))
      coords = 'geo'
    self.coords = coords

    # Set map projection limits and center point depending on hemisphere selection
    if lat_0 is None: 
      lat_0 = 90.
      if boundinglat: lat_0 = math.copysign(lat_0, boundinglat)
    if lon_0 is None: 
      lon_0 = -100.
      if self.coords == 'mag': 
        _, lon_0, _ = aacgm.aacgmConv(0., lon_0, 0., 0)
    if boundinglat:
      width = height = 2*111e3*( abs(lat_0 - boundinglat) )

    # Initialize map
    super(mapObjPlus, self).__init__(self, projection=projection, resolution=resolution, 
        lat_0=lat_0, lon_0=lon_0, width=width, height=height, **kwargs)

    # Add continents
    if coords is not 'mlt' or dateTime is not None:
      _ = self.drawcoastlines(linewidth=coastLineWidth)
      # self.drawmapboundary(fill_color=fillOceans)
      _ = self.fillcontinents(color=fillContinents, lake_color=fillLakes)

    # Add coordinate spec
    if showCoords:
      _ = text(self.urcrnrx, self.urcrnry, self._coordsDict[coords]+' coordinates', 
          rotation=-90., va='top', fontsize=8)

    # draw parallels and meridians.
    if grid:
      parallels = np.arange(-80.,81.,5.)
      out = self.drawparallels(parallels, color='.6', zorder=10)
      # label parallels on map
      if gridLabels: 
        lablon = int(self.llcrnrlon/10)*10
        rotate_label = lablon - lon_0 if lat_0 >= 0 else lon_0 - lablon + 180.
        x,y = basemap.Basemap.__call__(self, lablon*np.ones(parallels.shape), parallels)
        for ix,iy,ip in zip(x,y,parallels):
          if not self.xmin <= ix <= self.xmax: continue
          if not self.ymin <= iy <= self.ymax: continue
          _ = text(ix, iy, r"{:3.0f}$^\circ$".format(ip), 
              rotation=rotate_label, va='center', ha='center', zorder=10, color='.4')
      # label meridians on bottom and left
      meridians = np.arange(-180.,181.,20.)
      if gridLabels: 
        merLabels = [False,False,False,True]
      else: 
        merLabels = [False,False,False,False]
      # draw meridians
      out = self.drawmeridians(meridians, labels=merLabels, color='.6', zorder=10)

  
  def __call__(self, x, y, inverse=False, coords=None):
    from davitpy.models import aacgm
    from copy import deepcopy
    import numpy as np
    import inspect

    if coords is not None and coords not in self._coordsDict:
      print('Invalid coordinate system given in coords ({}): setting "{}"'.format(coords, self.coords))
      coords = None

    if coords and coords != self.coords:
      trans = coords+'-'+self.coords
      if trans in ['geo-mag','mag-geo']:
        flag = 0 if trans == 'geo-mag' else 1
        try:
          nx, ny = len(x), len(y)
          xt = np.array(x)
          yt = np.array(y)
          shape = xt.shape    
          y, x, _ = aacgm.aacgmConvArr(list(yt.flatten()), list(xt.flatten()), [0.]*nx, 2017, flag)
          x = np.array(x).reshape(shape)
          y = np.array(y).reshape(shape)
        except TypeError as e:
          y, x, _ = aacgm.aacgmConv(y, x, 0., 2017, flag)


    if self.coords is 'geo':
      return basemap.Basemap.__call__(self, x, y, inverse=inverse)

    elif self.coords is 'mag':
      try:
        callerFile, _, callerName = inspect.getouterframes(inspect.currentframe())[1][1:4]
      except: 
        return basemap.Basemap.__call__(self, x, y, inverse=inverse)
      if isinstance(y, float) and abs(y) == 90.:
        return basemap.Basemap.__call__(self, x, y, inverse=inverse)
      if 'mpl_toolkits' in callerFile and callerName is '_readboundarydata':
        if not inverse:
          try:
            nx, ny = len(x), len(y)
            x = np.array(x)
            y = np.array(y)
            shape = x.shape
            yout, xout, _ = aacgm.aacgmConvArr(list(y.flatten()), list(x.flatten()), [0.]*nx, 2017, 0)
            xout = np.array(xout).reshape(shape)
            yout = np.array(yout).reshape(shape)
          except TypeError:
            yout, xout, _ = aacgm.aacgmConv(y, x, 0., 0)
          return basemap.Basemap.__call__(self, xout, yout, inverse=inverse)
        else:
          return basemap.Basemap.__call__(self, x, y, inverse=inverse)
      else:
        return basemap.Basemap.__call__(self, x, y, inverse=inverse)

    elif self.coords is 'mlt':
      print('Not implemented')
      callerFile, _, callerName = inspect.getouterframes(inspect.currentframe())[1][1:4]


  def _readboundarydata(self, name, as_polygons=False):
    from davitpy.models import aacgm
    from copy import deepcopy
    import _geoslib
    import numpy as np

    if self.coords is 'mag':
      nPts = len(self._boundarypolyll.boundary[:, 0])
      lats, lons, _ = aacgm.aacgmConvArr(list(self._boundarypolyll.boundary[:, 1]), 
              list(self._boundarypolyll.boundary[:, 0]), 
              [0.]*nPts, 2017, 1)
      b = np.asarray([lons,lats]).T
      oldgeom = deepcopy(self._boundarypolyll)
      newgeom = _geoslib.Polygon(b).fix()
      self._boundarypolyll = newgeom
      out = basemap.Basemap._readboundarydata(self, name, as_polygons=as_polygons)
      self._boundarypolyll = oldgeom
      return out
    else: 
      return basemap.Basemap._readboundarydata(self, name, as_polygons=as_polygons)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -



def plot_all5(geographic=False, latspacing=10.):
    """
    Plot April 18th-22nd ephemeris tracks on same mapobj
    """
    fname_18th, idx_rev_18th, ellip_rev_18th = get_ottawa_data2("20160418")
    lons_18th, lats_18th, alts_18th, ephtimes_18th = get_rri_ephemeris(fname_18th)
    fname_19th, idx_rev_19th, ellip_rev_19th = get_ottawa_data2("20160419")
    lons_19th, lats_19th, alts_19th, ephtimes_19th = get_rri_ephemeris(fname_19th)
    fname_20th, idx_rev_20th, ellip_rev_20th = get_ottawa_data2("20160420")
    lons_20th, lats_20th, alts_20th, ephtimes_20th = get_rri_ephemeris(fname_20th)
    fname_21st, idx_rev_21st, ellip_rev_21st = get_ottawa_data2("20160421")
    lons_21st, lats_21st, alts_21st, ephtimes_21st = get_rri_ephemeris(fname_21st)
    fname_22nd, idx_rev_22nd, ellip_rev_22nd = get_ottawa_data2("20160422")
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
    elliplon_18th = float(lons_18th[ellip_rev_18th])
    elliplat_18th = float(lats_18th[ellip_rev_18th])
 
    shlon_19th = float(lons_19th[indx_shortest_19th])
    shlat_19th = float(lats_19th[indx_shortest_19th])
    invlon_19th = float(lons_19th[idx_rev_19th])
    invlat_19th = float(lats_19th[idx_rev_19th])
    elliplon_19th = float(lons_19th[ellip_rev_19th])
    elliplat_19th = float(lats_19th[ellip_rev_19th])

    shlon_20th = float(lons_20th[indx_shortest_20th])
    shlat_20th = float(lats_20th[indx_shortest_20th])
    invlon_20th = float(lons_20th[idx_rev_20th])
    invlat_20th = float(lats_20th[idx_rev_20th])
    elliplon_20th = float(lons_20th[ellip_rev_20th])
    elliplat_20th = float(lats_20th[ellip_rev_20th])

    shlon_21st = float(lons_21st[indx_shortest_21st])
    shlat_21st = float(lats_21st[indx_shortest_21st])
    invlon_21st = float(lons_21st[idx_rev_21st])
    invlat_21st = float(lats_21st[idx_rev_21st])
    elliplon_21st = float(lons_21st[ellip_rev_21st])
    elliplat_21st = float(lats_21st[ellip_rev_21st])

    shlon_22nd = float(lons_22nd[indx_shortest_22nd])
    shlat_22nd = float(lats_22nd[indx_shortest_22nd])
    invlon_22nd = float(lons_22nd[idx_rev_22nd])
    invlat_22nd = float(lats_22nd[idx_rev_22nd])
    elliplon_22nd = float(lons_22nd[ellip_rev_22nd])
    elliplat_22nd = float(lats_22nd[ellip_rev_22nd])
   
    # A different font for the legend etc. might be nice
    fig = plt.figure(1, figsize=(14,8))
    font = {'fontname':'Computer Modern'}

    # Option 1: use my custom mapObjPlus class
    #m = mapObjPlus(lat_0=57.0, lon_0=5.0, width=18e3*180, height=25e3*90, coords='mag',resolution='i',datetime=times_20th[0])
  
    # Option 2: use a modified version of the classic thing (modify plotUtils.mapObj locally to have divisions every 5 deg 
    if geographic==True:
        m = plotUtils.mapObj(lat_0=46.0, lon_0=-75.0, width=18e3*180, height=25e3*90, coords='geo',resolution='i',datetime=times_20th[0], gridLatRes=latspacing, fillOceans=(1.,1.,1))
    else:
        m = plotUtils.mapObj(lat_0=57.0, lon_0=5.0, width=18e3*180, height=25e3*90, coords='mag',resolution='i',datetime=times_20th[0], gridLatRes=latspacing, fillOceans=(1.,1.,1))

    # Option 3: use the regular classic davitpy.plotUtils.mapObj and try to redraw more parallels
    #m = plotUtils.mapObj(lat_0=57.0, lon_0=5.0, width=18e3*180, height=25e3*90, coords='mag',resolution='i',datetime=times_20th[0])
    #parallels = np.arange(0., 81., 5.)
    #m.drawparallels(parallels)

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

if __name__=="__main__": 
    plot_all5(geographic=True, latspacing=10.)
    plot_all5(geographic=True, latspacing=5.)
