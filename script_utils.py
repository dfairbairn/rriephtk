"""
File: 'script-utils.py'
Description: Contains utilities used in the RRI Conjunction Finder script.
    These include functions to more conveniently grab FOV data from Davitpy,
    plotting functions that make use of Matplotlib, etc.
Author: David Fairbairn
Date: June 2016


"""
import matplotlib.pyplot as plt
import subprocess
import numpy as np
import datetime as dt #A wrapper class for 'file' which tracks line numbers. 

class FileLineWrapper(object):
    """
    A wrapper class for 'file' which tracks line numbers.

    **ATTRIBUTES**
        f (file object): the file object being wrapped.
        line (int): tracks current line number
        byte_offs (int):  
    """
    def __init__(self, f):
        """ Constructor for filewrapper object. """
        self.f = f
        self.line = 0
        
        # To allow skipping directly to lines
        self.line_offs = []
        offset = 0
        self.f.seek(0)
        for line in f:
            #print line
            self.line_offs.append(offset)
            offset += line.__len__()
        self.f.seek(0)
    def close(self):
        """ Closer function for filewrapper object. """
        return self.f.close()
    def readline(self):
        """  """
        self.line += 1
        return self.f.readline()
    def seekline(self,line_num):
        """  """ 
        self.f.seek(self.line_offs[line_num - 1])

def ephems_to_datetime(ephem_times):
    """
    This function allows a whole array of ephemeris times to be conveniently
    converted in the same manner as with the ephem_to_datetime() function

    ** ARGS **
        ephem_times (numpy.array): array of Ephemeris MET (truncated JD times)

    ** RETURNS **
        times (list): list of datetime objects
    """
    assert isinstance(ephem_times, np.ndarray), "Not an array"
 
    # i) Check the seconds between May 24 1968 (ephem MET) and Jan 1 1970 (neg number)
    t_off = subprocess.check_output(["date", "--date=1968-05-24 0:00:00", "+%s", "-u"])
    t_off = float(t_off.split("\n",1)[0]) # extract the integer value

    # ii) Do the math.
    ephtimes = ephem_times + t_off
    times = []
    for i in range(np.size(ephtimes)):
        times.append(dt.datetime.utcfromtimestamp(ephtimes[i]))
    return times

def ephem_to_datetime(ephem):
    """
    This function exists in order to conveniently convert the weird ephemeris
    time data format from the EPOP RRI instrument into a datetime object.

    *note*: makes use of the bash command 'date', so make sure you have it.

    ** ARGS **
        ephem (int): a large integer showing # of seconds since May 24 1968
                    (must be greater than 0, smaller than sys.maxint)

    ** RETURNS **
        dtime (datetime.datetime): a datetime object for the time given in ephem 
    """    
    assert ephem >= 0
    assert isinstance(ephem, int)

    # i) Check # seconds between May 24 1968 (ephem MET) and Jan 1 1970 (neg number)
    # -u parameter required to ensure offset is in UTC. +%s specifies output in seconds.
    t_off = subprocess.check_output(["date", "--date=1968-05-24 0:00:00", "+%s", "-u"])
    t_off = float(t_off.split("\n",1)[0]) # extract the integer value
  
    # ii) Do the math. 
    dtime = dt.datetime.utcfromtimestamp(t_off + float(ephem))    
    return dtime


#TODO: plot the back lobe?
def plot_fov_sat(fovname, date, ephem_lons, ephem_lats, suppress_show=False):
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
    fov = get_fov_by_name(fovname)
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

    plt.title(fovname + " Radar FOV vs RRI Ephemeris for " + date.__str__())
    plt.xlabel('Geographic Longitude (degrees)')
    plt.ylabel('Geographic Latitude (degrees)')

    # So that the formats match, I will ensure months and days are padded for
    # the output figure's name as well.
    month = "0" + str(date.month) if str(date.month).__len__() == 1 else str(date.month)
    day = "0" + str(date.day) if str(date.day).__len__() == 1 else str(date.day)
    hr = "0" + str(date.hour) if str(date.hour).__len__() == 1 else str(date.hour)
    mn = "0" + str(date.minute) if str(date.minute).__len__() == 1 else str(date.minute)

    plt.plot(ephemlons,ephemlats,'r',label="RRI Ephemeris")
    plt.legend()
    plt.savefig("./data/output/"+str(date.year)+month+day+"_"+hr+"h"+mn+"_"+fovname) 
    if suppress_show==False:
        plt.show()
    
def get_fov_by_name(name):
    """
    This function shortcuts the FOV creation process by handling assumptions I
    typically make when creating FOV objects.

    ** ARGS **
        name (string): the name (as DaVitPy understands them) of a SuperDARN 
            radar whose FOV is to be returned.
    
    ** RETURNS **
        fov (pydarn.radar.radFov.fov): a (front) FOV object for the SuperDARN
            site with the given name.
    """
    from davitpy import pydarn
    nw = pydarn.radar.network()
    rad = nw.getRadarByName(name)
    assert (rad != False) #if name is illegal, getRadarByName returns False
    site = pydarn.radar.site(radId = rad.id)
    fov = pydarn.radar.radFov.fov(site=site,altitude=300.0,model='IS',coords='geo',ngates=75)
    return fov

def two_pad(in_time):
    """ 
    Takes in a number of 1 or 2 digits, returns a string of two digits. 
    """
    assert isinstance(in_time,int)
    assert (in_time < 100 and in_time >= 0)
    return "0" + str(in_time) if str(in_time).__len__() == 1 else str(in_time)


"""

TESTING

"""
if __name__ == "__main__":

    # Note: Can't test illegal inputs because I'm not using exceptions,
    # but the assertions in the methods should hopefully be fine for now.

    # Testing FileLineWrapper
    f = open('./script_utils.py','r')
    f.readline()
    f.readline()
    offs = f.tell()
    line = f.readline()
    fw = FileLineWrapper(f)
    assert(fw.line_offs.__len__() > 0)
    assert(fw.line_offs[2] == offs)
    fw.seekline(3) 
    assert(fw.readline() == line)
    linen1 = fw.line
    fw.readline()
    assert(linen1 + 1 == fw.line)

    # Testing ephem_to_datetime():
    tst_dt1 = dt.datetime(1968,5,24)
    assert(ephem_to_datetime(0) == tst_dt1)
    tst_dt2 = dt.datetime(2014,7,8,1,15,9)
    num2 = 1455498909
    assert(ephem_to_datetime(num2) == tst_dt2)
 
    # Testing ephems_to_datetime():
    assert(ephems_to_datetime(np.array([0,num2])) == [tst_dt1, tst_dt2])

    # Testing two_pad():
    assert(two_pad(3) == "03")
    assert(two_pad(13) == "13")

    # Testing get_fov_by_name():
    from davitpy import pydarn
    assert isinstance(get_fov_by_name("Saskatoon"),pydarn.radar.radFov.fov)

    # Testing plot_fov_sat():
    # Check for an output file?
    # We can do this some day in the future...
