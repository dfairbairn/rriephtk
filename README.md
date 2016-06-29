# RRI-Conjunction-Script
A tool which locates relevant SuperDARN signals which might have reached the RRI instrument.

This is a work-in-progress.


**Install:**
This tool has *not* been made particularly portable yet. An automatic install
process does not exist: you'll need to download and install DaVitPy separately, 
which in turn requires Numpy and Matplotlib which are used liberally in this 
script.

The script will create several subdirectories for mounting and saving data:
./data, ./data/output, ./data/remote. It initially checks to ensure they are
present, and if not, creates them.

In addition, some of the data is retrieved from a server from which 
authentication is done by the user, so unless you have access to it
already, that won't work either. 

**Code:**
The important files are 'script.py', 'script-utils.py', and 'range_cells.py'.
The files 'radar_finder.py' and 'rgCoords.py' contain different options for 
radar intersection determination which were used or considered previously, but
are not currently in use. They will likely be removed soon.

**Usage:**
If you happen to have a computer with all the correct software and you have
access to the data server, the script can be used by calling it at the 
command line with the desired RRI data file as an additional argument:

If pwd is script directory (or optionally use the -i option to tell Python
to enter interactive mode after running the script):
> python script.py ./data/RRI_20160401_072744_033718_lv1_v2.h5

> python -i script.py ./data/RRI_20160401_072744_033718_lv1_v2.h5
