""" file: magnet_data.py
author: David Fairbairn
date: September 2016

Looking at Magnetometer data from ePOP during the Ottawa passes of April 2016.

"""
from data_utils import *
from ottawa_plots import *
import h5py

sample_mgf_datafile = "data/mgf/MGF_20160418_222505_224033_V_01_00_00.1sps.GEI.lv3"

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

def satframe2ned(sc_vecs,ram_dirs,yaws,pitchs,rolls):
    """
    This function takes coordinates expressed in terms of directions with 
    respect to the satellite body. The x,y,z coordinates would normally 
    correspond to regular spacecraft coordinates, but in this case the body
    of the spacecraft is assumed to be oriented away from its ram direction, 
    which is corrected for.

    This correction is achieved by accounting for the craft's yaw, pitch, roll
    values at each ephemeris point.

    Note: yaw, pitch, roll are angles of rotation around three 'SC' coordinate
    axes. 'X' is in the ram direction, 'Z' is the axis from the craft body to the
    Earth's centre, and 'Y' is 'Z' cross 'X'.

    The assumption in this function is that the coordinates in sc_vecs are 
    not actually yet in proper SC/('X','Y','Z') components yet - that they are
    tilted around with the satellite body.

    *** PARAMS ***
    sc_vecs: the coords in terms of the spacecraft body (x=ram, z=nadir)
    ram_dirs: the ram direction in N-E-D components
    yaws: rotation around spacecraft Z axis in degrees
    pitchs: rotation around spacecraft Y axis in degrees
    rolls: rotation around spacecraft X axis in degrees

    *** RETURNS ***
    outs: output vector (should be in same units as input)
 
    **!!Standards Dilemma!!**
    Word of Gareth is that x is ram dir, z is nadir dir, and y is Z cross X.
    But then if satellite has vertical component (which it does, though its 
    small), this isn't a right handed coord system. Possible solutions are:
    - ignore it and hope the error is small (should be)
    - define x to just be ram direction in North and East directions (not down)

    Initially, I proceed with solution approach #1 for simplicity 
 
    """
    # Step 1. Compute Spacecraft X Y Z directions from the ram direction
    xdirs = [ram_dirs[i]/np.linalg.norm(ram_dirs[i]) for i in range(len(ram_dirs))]

    print "First normalized ram dir (x dir): ",xdirs[0] #reminder to self: check if its already normalized
    zdirs = np.array([(0.,0.,1.) for i in range(xdirs.__len__())]) #Just down
    ydirs = np.cross(zdirs,xdirs)



    # Step 2. Reverse rotate by amounts described in roll, pitch, yaw
    outs = []
    for i in range(xdirs.__len__()):
        roll_rot = rotation_matrix(xdirs[i],np.deg2rad(-rolls[i]))
        pitch_rot = rotation_matrix(ydirs[i], np.deg2rad(-pitchs[i]))
        yaw_rot = rotation_matrix(zdirs[i], np.deg2rad(-yaws[i]))
        #print "Roll rotation matrix magnitude: ",np.linalg.norm(roll_rot) 
        #print "Pitch rotation matrix magnitude: ",np.linalg.norm(pitch_rot) 
        #print "Yaw rotation matrix magnitude: ",np.linalg.norm(yaw_rot) 
        """
        intermed1 = np.dot(roll_rot, sc_vecs[i])
        intermed2 = np.dot(pitch_rot, intermed1)
        intermed3 = np.dot(yaw_rot, intermed2)
        """
        intermed3 = sc_vecs[i] # Temporarily looking at what happens if I don't bother accounting for yaw/roll/pitch
        A = np.array((xdirs[i],ydirs[i],zdirs[i]))
        Ainv = np.linalg.inv(A)

        #print "Spatial conversion matrix: ",A
        #print "Spatial conversion matrix magnitude: ",np.linalg.norm(A)
        #exit_rri()
        out = np.dot(intermed3, Ainv)
        outs.append(out)

    return np.array(outs)

def sc2ned(sc_vecs,ram_dirs):
    """
    Takes some space-craft based coordinates and the spacecraft's ram direction
    in North-East-Down components and converts another set of vectors written
    in terms of spacecraft coordinates, and converts those to NED.

    *** PARAMS ***
    sc_vecs: the coords in terms of the spacecraft frame (x=ram, z=nadir)
    ram_dirs: the ram direction in N-E-D components

    *** RETURNS ***
    outs: output vector (should be in same units as input)
 
    **!!Standards Dilemma!!**
    Word of Gareth is that x is ram dir, z is nadir dir, and y is Z cross X.
    But then if satellite has vertical component (which it does, though its 
    small), this isn't a right handed coord system. Possible solutions are:
    - ignore it and hope the error is small (should be)
    - define x to just be ram direction in North and East directions (not down)

    Initially, I proceed with solution approach #1 for simplicity 
 
    """
    # Step 1. Compute Spacecraft X Y Z directions from the ram direction
    xdirs = [ram_dirs[i]/np.linalg.norm(ram_dirs[i]) for i in range(len(ram_dirs))]

    print "First normalized ram dir (x dir): ",xdirs[0] #reminder to self: check if its already normalized
    zdirs = np.array([(0.,0.,1.) for i in range(xdirs.__len__())]) #Just down
    ydirs = np.cross(zdirs,xdirs)


    outs = []
    for i in range(xdirs.__len__()):

        intermed = sc_vecs[i] # Temporarily looking at what happens if I don't bother accounting for yaw/roll/pitch
        A = np.array((xdirs[i],ydirs[i],zdirs[i]))
        Ainv = np.linalg.inv(A)

        #print "Spatial conversion matrix: ",A
        #print "Spatial conversion matrix magnitude: ",np.linalg.norm(A)
        #exit_rri()
        out = np.dot(intermed, Ainv)
        outs.append(out)
    return np.array(outs)

def ned2sc(ned_vecs,ram_dirs):
    """
    Also need a way to compute SC coords from NED coords. This function is the
    inverse of sc2ned3
    """
    # Step 1. Compute Spacecraft X Y Z directions from the ram direction
    xdirs = [ram_dirs[i]/np.linalg.norm(ram_dirs[i]) for i in range(len(ram_dirs))]

    print "First normalized ram dir (x dir): ",xdirs[0] #reminder to self: check if its already normalized
    zdirs = np.array([(0.,0.,1.) for i in range(xdirs.__len__())]) #Just down
    ydirs = np.cross(zdirs,xdirs)

    outs = []
    for i in range(xdirs.__len__()):

        intermed = ned_vecs[i] # Temporarily looking at what happens if I don't bother accounting for yaw/roll/pitch
        A = np.array((xdirs[i],ydirs[i],zdirs[i]))

        #print "Spatial conversion matrix: ",A
        #print "Spatial conversion matrix magnitude: ",np.linalg.norm(A)
        #exit_rri()
        out = np.dot(A,intermed)
        outs.append(out)
    return np.array(outs)
   

def cmp_igrf_magnetometer(fname=sample_mgf_datafile, date_string="20160418"):
    """
    I use this to compare the IGRF model to the recorded magnetic field by the
    MGF instrument onboard ePOP.
    """
    # Set this to correspond to the mgf file at the top until mgf file selection is possible
    
    datpath,datname = initialize_data()
    rrifname,index_reversal = get_ottawa_data(date_string)
    
    # Get RRI ephemeris data together so that we can remove effect of spacecraft direction
    lons,lats,alts,ephtimes,mlons,mlats,mlts,pitchs,yaws,rolls = get_rri_ephemeris_full(rrifname)
    ephtimes = np.array([ round(e) for e in ephtimes]) # crucial for comparing mgf and rri times
    vs,dists = get_ramdirs(lons, lats, alts, ephtimes)

    # Calculate IGRF at each ephemeris point
    # import os
    # import sys
    # Redirect stdout so that IGRF can't print to the screen and spam me
    # sys.stdout = open(os.devnull, "w") # turns out it doesnt work because IGRF stuff is in fortran/handled separately
    B_igrf,kvecs,angles = get_kb_ottawa_angle(lons,lats,alts,ephtimes) 
    # sys.stdout = sys.__stdout__ # Reconnect stdout

    # TMP change! Leaving both arrays of B field measurements unchanged and playing with spacepy on them
    #B_igrf = ned2sc(np.array(B_igrf),vs)  #Converting from NED to SC!
    B_igrf = np.array(B_igrf)
 
    # Acquire the MGF data
    bscx, bscy, bscz, ephtimes_bsc = read_mgf_file(fname)
    B_mgf_intermediate = [ (bscx[i],bscy[i],bscz[i]) for i in range(len(bscx))]  
    print "Intermediate B_mgf first entry:\n",B_mgf_intermediate[0]
    print "Magnitude: ",np.linalg.norm(B_mgf_intermediate[0])

    # Need to compare the time of the different data sets 
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

    B_mgf = np.array(B_mgf_intermediate)

    return B_mgf,B_igrf

def plot_comparison(fname=sample_mgf_datafile, date_string="20160418"):
    """
    Function to plot the components of MGF and IGRF B fields against each other.
    Note: Currently the IGRF might not be correctly converted so the components
    won't be similar to each other at all.
    """
    rri_fname,index_reversal = get_ottawa_data(date_string)
    lons,lats,alts,ephtimes,mlons,mlats,mlts,pitchs,yaws,rolls = get_rri_ephemeris_full(rri_fname)
   
    B_mgf,B_igrf = cmp_igrf_magnetometer(fname,date_string)
    plt.plot(B_mgf[:,0],'b',label="MGF SC_X Component")
    plt.plot(B_igrf[:,0],'r',label="IGRF SC_X Component")
    plt.legend()
    plt.title("Comparison of B Component in SC_X Direction")
    ephem_ticks(lons,lats,alts,ephtimes,mlons,mlats,mlts)
    plt.show()

    plt.plot(B_mgf[:,1],'b',label="MGF SC_Y Component")
    plt.plot(B_igrf[:,1],'r',label="IGRF SC_Y Component")
    plt.legend()
    plt.title("Comparison of B Component in SC_Y Direction")
    ephem_ticks(lons,lats,alts,ephtimes,mlons,mlats,mlts)
    plt.show()

    plt.plot(B_mgf[:,2],'b',label="MGF SC_Z Component")
    plt.plot(B_igrf[:,2],'r',label="IGRF SC_Z Component")
    plt.legend()
    plt.title("Comparison of B Component in SC_Z Direction")
    ephem_ticks(lons,lats,alts,ephtimes,mlons,mlats,mlts)
    plt.show()

if __name__ == "__main__":
    
    B_mgf,B_igrf = cmp_igrf_magnetometer()
    print "MGF B field magnitude: ",np.linalg.norm(B_mgf[0]),B_mgf[0]
    print "IGRF B field magnitude: ",np.linalg.norm(B_igrf[0]),B_igrf[0]
    #plot_comparison()

    """
    bscx, bscy, bscz, ephtimes_bsc = read_mgf_file(sample_mgf_datafile)
    B_mgf = [ (bscx[i],bscy[i],bscz[i]) for i in range(len(bscx))]  
    #print B_mgf
 
    date_string = "20160418"
    fname, index_reversal = get_ottawa_data(date_string)
    lons,lats,alts,ephtimes,mlons,mlats,mlts,pitchs,yaws,rolls = get_rri_ephemeris_full(fname)
    vs,dists = get_ramdirs(lons,lats,alts,ephtimes)   
    
    print B_mgf[0] 
    out = sc2ned2(B_mgf,vs,yaws,pitchs,rolls)
    """
