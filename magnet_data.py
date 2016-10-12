""" file: magnet_data.py
author: David Fairbairn
date: September 2016

Looking at Magnetometer data from ePOP during the Ottawa passes of April 2016.

"""
from data_utils import *
from ottawa_plots import *
import h5py

sample_mgf_datafile = "data/mgf/MGF_20160418_222505_224033_V_01_00_00.1sps.SC.lv3"

def read_mgf_file(fname):
    # if doesn't exist, exit with error
    
    # else try to read in the columns as we expect them

    #f = FileLineWrapper(open(fname,'r'))
    # FileLineWrapper functionality (seeking/returning line #'s) isn't really
    # necessary for reading in file data once - so I don't bother using it
    ephtimes = []
    Bscx = []
    Bscy = []
    Bscz = []

    f = open(fname,'r')
    for ln in f:
        spl = ln.split()
        
        if spl[0]=='LV3_DESCRPT':
            continue 
        ## TODO: Make this check work for 5 decimal doubles (alt to isdigit()?)
        #if not spl[1].isdigit():
        #    continue

        # No need to retrieve other parts of data
        ephtimes.append(float(spl[1])) # Ephem Times
        Bscx.append(float(spl[2])) # B_SCx
        Bscy.append(float(spl[3])) # B_SCy
        Bscz.append(float(spl[4])) # B_SCz
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

    #A = np.array([scx,scy,scz])
    A = np.array((x,y,z))
    #print "A:\n",A
    Ainv = np.linalg.inv(A)

    #print "scx:\n",scx
    #print "scy:\n",scy
    #print "scz:\n",scz
    #print "Ainv:\n",Ainv

    (outN,outE,outD) = np.dot((float(scx),float(scy),float(scz)),Ainv)
    return (outN,outE,outD)

def cmp_igrf_magnetometer():

    # Set this to correspond to the mgf file at the top until mgf file selection is possible
    date_string = "20160418"
    datpath,datname = initialize_data()
    fname,index_reversal = get_ottawa_data(date_string)
    
    # Get ephemeris data together so that 
    lons,lats,alts,ephtimes = get_rri_ephemeris(fname)
    ephtimes = np.array([ round(e) for e in ephtimes]) # crucial for comparing mgf and rri times
    v,dists = get_ramdirs(lons, lats, alts, ephtimes)
    bscx, bscy, bscz, ephtimes_bsc = read_mgf_file(sample_mgf_datafile)

    times_rri =  ephems_to_datetime(ephtimes)
    times_mgf = ephems_to_datetime(np.array(ephtimes_bsc))

    print "times_mgf's first few entries look like:\n",times_mgf[0:3]
    print "times_rri's first few entries look like:\n",times_rri[0:3] 
    print "Length of mgf times:\n",len(ephtimes_bsc)
    print "Length of rri times:\n",len(ephtimes)
    try:
        i_rristart = times_mgf.index(times_rri[0]) 
        print "Index of mgf data at which point rri starts taking data:\t",i_rristart
        print "times_mgf_iso[i_rristart]:\t",times_mgf[i_rristart]
        print "times_rri_iso[0]:\t",times_rri[0]
    except ValueError:
        print "Failed to find where RRI starts in MGF data." 
        print times_rri_iso

    #TODO: change sc2ned so that we don't need this expensive loop    
    bscN = []
    bscE = []
    bscD = []
    for i in range(len(v)):
        vp = v[i]
        print str(i)+"th MGF times entry:\t",times_mgf[i+i_rristart]
        print str(i)+"th RRI times entry:\t",times_rri[i]
        bscn_tmp, bsce_tmp, bscd_tmp = sc2ned(bscx[i+i_rristart],bscy[i+i_rristart],bscz[i+i_rristart],vp[0],vp[1],vp[2])
        bscN.append(bscn_tmp)
        bscE.append(bsce_tmp)
        bscD.append(bscd_tmp)

    B_mgf = np.array([ (bscN[i],bscE[i],bscD[i]) for i in range(len(bscN))])
    B_igrf,kvecs,angles = get_kb_ottawa_angle(lons,lats,alts,ephtimes)
    
    return B_mgf,B_igrf

if __name__ == "__main__":
    B_mgf,B_igrf = cmp_igrf_magnetometer()
    
    B_N_mgf = [b[0] for b in B_mgf]
    B_E_mgf = [b[1] for b in B_mgf]
    B_D_mgf = [b[2] for b in B_mgf]
    B_N_igrf = [b[0] for b in B_igrf]
    B_E_igrf = [b[1] for b in B_igrf]
    B_D_igrf = [b[2] for b in B_igrf] 


