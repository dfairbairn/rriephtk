""" file: magnet_data.py
author: David Fairbairn
date: September 2016

Looking at Magnetometer data from ePOP during the Ottawa passes of April 2016.

"""
from data_utils import *
from ottawa_plots import *
import h5py

sample_mgf_datafile = "data/mgf/MGF_20160420_055205_061332_V_01_00_00.1sps.SC.lv3"

def read_mgf_file(fname):
    # if doesn't exist, exit with error
    
    # else try to read in the columns as we expect them

    #f = FileLineWrapper(open(fname,'r'))
    # FileLineWrapper functionality (seeking/returning line #'s) isn't really
    # necessary for reading in file data once - so I don't bother using it
    ephtimes = Bscx = Bscy = Bscz = []

    f = open(fname,'r')
    for ln in f:
        spl = ln.split()
        #print spl
        
        if not spl[1].isnumeric():
            continue

        # No need to retrieve other parts of data
        ephtimes.append(spl[1]) # Ephem Times
        Bscx.append(spl[2]) # B_SCx
        Bscy.append(spl[3]) # B_SCy
        Bscz.append(spl[4]) # B_SCz
    return (Bscx, Bscy, Bscz, ephtimes) 
    
def sc2ned(scx,scy,scz,ramN,ramE,ramD):
    """
    Takes some space-craft based coordinates and the spacecraft's ram direction
    in North-East-Down components, and converts the space-craft coords to NED.

    *** PARAMS ***
    scx, scy, scz: the coords in terms of the spacecraft frame (x=ram, z=nadir)
    ramN,ramE,ramD: the ram direction in N-E-D components

    *** RETURNS ***
    outN, outE, outD: converts
 
    **!!Standards Dilemma!!**
    Word of Gareth is that x is ram dir, z is nadir dir, and y is Z cross X.
    But then if satellite has vertical component (which it does, though its 
    small), this isn't a right handed coord system. Possible solutions are:
    - ignore it and hope the error is small (should be)
    - define x to just be ram direction in North and East directions (not down)

    InitiallyFor , I proceed with solution approach #1 for simplicity 
   
    """
    
    #TODO: validate inputs  
    x = (ramN, ramE, ramD)
    #print scx #reminder to self: check if its already normalized
    z = (0.,0.,1.) #Just down
    y = np.cross(z,x)

    A = np.array([scx,scy,scz])
    Ainv = np.linalg.inv(A)
    outN = np.dot(scx,Ainv)
    outE = np.dot(scy,Ainv)
    outD = np.dot(scz,Ainv) 
    return (outN,outE,outD)


def cmp_igrf_magnetometer():
    dat1,dat2 = initialize_data()
    f = h5py.File(dat2)
    return 

if __name__ == "__main__":
    # Set this to correspond to the mgf file at the top until mgf file selection is possible
    date_string = "20160420"
    datpath,datname = initialize_data()
    fname,index_reversal = get_ottawa_data(date_string)
    
    # Get ephemeris data together so that 
    #
    lons,lats,alts,ephtimes = get_rri_ephemeris(fname)
    v,dists = get_ramdirs(lons, lats, alts, ephtimes)
    bscx, bscy, bscz, ephtimes_bsc = read_mgf_file(sample_mgf_datafile)

    print len(ephtimes_bsc)
    print len(ephtimes)
    print "ephtimes_bsc's first few entries look like:\n",ephtimes_bsc[0:3]
    print "ephtimes's first few entries look like:\n",ephtimes[0:3] 
    
    print ephtimes_bsc[1].isnumeric()


    #TODO: chagne sc2ned so that we don't need this expensive loop    
    bscN = bscE = bscD = []
    for i in range(len(bscx)):
        bscn_tmp, bsce_tmp, bscd_tmp = sc2ned(bscx,bscy,bscz,v[0],v[1],v[2])
        bscN.append(bscn_tmp)
        bscE.append(bsce_tmp)
        bscD.append(bscd_tmp)
