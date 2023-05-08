# RRI-EphTK
The RRI Ephemeris Toolkit (RRI-EphTK) is a mini-toolkit for performing various analyses 
on .hdf5 files from the Radio Receiver Instrument (RRI) in the ePOP scientific payload of 
the CASSIOPE satellite.

This started out as a script looking at *conjunctions* of the satellite passes coinciding 
with passing through the fields of view of the SuperDARN radars (see conjunctions). Now 
it includes plotting tools, analytical tools for estimating various geometric features 
from transmitter to RRI, and some ionospheric parameter modeling.


## Installation
This tool has *not* been made particularly portable yet. A surefire automatic install
process does not exist. However, you _ought_ to be able to get it going by:

  i) Cloning this repository 
  > \> git clone http://github.com/dfairbairn/rriephtk

  ii) Setting up a virtual environment ( http://docs.python-guide.org/en/latest/dev/virtualenvs/ )
  > \> cd rriephtk
  > \> virtualenv env
  > \> source env/bin/activate

  iii) and then installing from the requirements.txt document here.
  > \> pip install -r requirements.txt

The script will create several subdirectories for mounting and saving data:
./data, ./data/output, ./data/remote. It initially (when data_utils.initialize_data 
is called) checks to ensure they are present, and if not, creates them.

In addition, some of the data is retrieved from a server from which 
authentication is done by the user, so unless you have access to it
already, that won't work either. 

## Usage

_Conjunctions_
If you happen to have a computer with all the correct software and you have
access to the data server, the script can be used by calling it at the 
command line with the desired RRI data file as an additional argument:

If pwd is rriephtk top directory (or optionally use the -i option to tell Python
to enter interactive mode after running the script):
> python rriephtk/conjunctions/conjunctions.py ./data/RRI_20160401_072744_033718_lv1_v2.h5

> python -i rriephtk/conjunctions/conjunctions.py ./data/RRI_20160401_072744_033718_lv1_v2.h5

Reading the documentation and using the functions directly could work too, since it's not super externally usable yet!


## Structure

The 'analysis' package contains 'magnet_data' and 'analysis_tools' for performing
estimation, simulation, and analysis.

A package for plots ('plots') will allow a few configurations of plots of 
CASSIOPE ground tracks (including some with SuperDARN FOV's overlaid). Various
parameters that are estimated in the 'analysis' package can be plotted using these
functions.

The 'conjunctions.conjunctions' package contains data which runs through coordinates 
of RRI geographic positions, determines which SuperDARN radars may be relevant,
and summarizes them and the contents of their 'errlog' information files from
that time period. Further study is available using the 'time_align' module.

### Analysis Tools
This package is intended to function purely as a library that can be imported.
Directly running 'analysis_tools' as above will result in a variety of automatic tests and 
demonstrations to be run. 

### Plotting_
Library functionality is also intended for this package, with the exception of the 'fig1_plots.py' file,
which can be run on its own (if you have the specified RRI data files) to plot 5 ground tracks near Ottawa
from April 2016.

### Utils
Importing and using the data_utils package is necessary for use of this package. Use initialize_data()
to ensure a data directory structure is set up and the script is connected to the remote server containing
various data sources.
