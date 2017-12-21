"""
file: 'check_rri_conditions.py'
description:
    Really analyzes the conditions from SuperDARN's POV surrounding an RRI pass.
author: David Fairbairn
date: August 15th 2017

"""
import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
import datetime as dt

import rriephtk_path
import rriephtk.utils.data_utils as datautils
import rriephtk.plotting as plotting
import rriephtk.conjunctions.conjunctions as conj

def radio_data(fname): 
    import h5py 
    f = h5py.File(fname, 'r')
    rridat = f['RRI Data']
    rm1 = rridat['Radio Data Monopole 1 (mV)'].value
    rm2 = rridat['Radio Data Monopole 2 (mV)'].value
    rm3 = rridat['Radio Data Monopole 3 (mV)'].value
    rm4 = rridat['Radio Data Monopole 4 (mV)'].value

    n_groups = len(rm1)
    # Finding and removing the nan's   
    for i in range(n_groups):
        datautils.update_progress(float(i)/n_groups)
        group = rm1[i]
        if all([ not np.isnan(g) for g in group]):
            pass
        else:
            print("\nRM1[{0}]: {1}".format(i, group))
            print("RM2[{0}]: {1}".format(i, rm2[i]))
            print("RM3[{0}]: {1}".format(i, rm3[i]))
            print("RM4[{0}]: {1}".format(i, rm4[i]))
            for j in range(len(group)):
                rm1[i][j] = 0.
                rm2[i][j] = 0.
                rm3[i][j] = 0.
                rm4[i][j] = 0.

    return rm1, rm2, rm3, rm4
    
    # Using every 29-length bundle to calculate stokes parameters
    #chis = np.array(n_groups)
    #psis = np.array(n_groups)
    chis = []
    psis = []
    for i in range(n_groups):
        datautils.update_progress(float(i)/n_groups)
        rm1group = rm1[i]
        rm2group = rm2[i]
        rm3group = rm3[i]
        rm4group = rm4[i]

        # TODO: Vectorize
        phi1 = np.arctan2(rm2group, rm1group)
        phi3 = np.arctan2(rm4group, rm3group)
        delta_phi = phi1 - phi3
        I = np.mean(np.power(rm1group, 2) ) + np.mean(np.power(rm3group, 2))
        Q = np.mean(np.power(rm1group, 2) ) - np.mean(np.power(rm3group, 2))
        U = 2*np.mean(rm1group*rm3group*np.cos(delta_phi))
        V = 2*np.mean(rm1group*rm3group*np.sin(delta_phi))

        Ip = np.sqrt(Q**2 + U**2 + V**2)
        m  = Ip/I
        psi = np.arctan2(U, Q)/2.
        psi_deg = np.rad2deg(psi)
        psis.append(psi)
        chi = np.arcsin(np.array([V]), np.array([Ip]))[0]/2.
        chi_deg = np.rad2deg(psi)
        chis.append(chi)

    psis = np.array(psis)
    chis = np.array(chis)
    return psis, chis

def create_pulseseqs(freq=62500.33933):
    """
    Makes the 7- and 8-pulse sequences that happen in SuperDARN

    [:param freq:]: sampling frequency that would be detecting the sequence
    
    """
    pulse_diffs_8 = [0.0210, 0.012, 0.003, 0.0045, 0.006, 0.0165, 0.0015]
    pulse_diffs_7 = [0.0216, 0.0072, 0.0192, 0.0048, 0.0096, 0.0024]
    pseq_length = 0.07

    n_samples = pseq_length*freq
    pulse_t = np.arange(0, pseq_length, pseq_length/n_samples)
    pulse_7 = 0.*pulse_t
    pulse_8 = 0.*pulse_t 
  
    pulse_7[0] = 1. 
    pulse_8[0] = 1.
    for i, e in enumerate(pulse_diffs_7): 
        offset_time = sum(pulse_diffs_7[:i+1])
        print("Pulse tdiff: {0}, Overall time into sequence: {1}".format(e, offset_time))
        print("Index to insert at: {0}".format(int(offset_time*freq)))
        pulse_7[int(offset_time*freq)] = 1
         
    for i, e in enumerate(pulse_diffs_8):
        offset_time = sum(pulse_diffs_8[:i+1])
        pulse_8[int(offset_time*freq)] = 1

    return pulse_t, pulse_7, pulse_8 
    



if __name__=="__main__":
    print("Enter the RRI filename to use")
    fname_rri = raw_input() 
    if not os.path.isfile(fname_rri):
        print("Invalid file, proceeding with default 20170419...")
        fname_rri = './data/RRI_20170419_225444_230041_lv1_v3.h5'
    lons, lats, alts, ephtimes = datautils.get_rri_ephemeris(fname_rri)
    times = datautils.ephems_to_datetime(ephtimes)
    t0 = times[0]
    ''' 
    #radar_conjs = conj.get_conjunctions(fname_rri)
    codes = ['fhe']
    if 'kap' in codes:
        codes.remove('kap')
 
    for code in codes:
        for i in np.arange(-10, 15):
            tdelt = dt.timedelta(0, 0, 0, 0, i)
            t = t0 + tdelt
            print(t)

            plotting.plot_fan_rri(lons, lats, alts, t, [code], save=False)
            break
    '''
