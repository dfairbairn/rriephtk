# RRI-Conjunction-Script
A tool which locates relevant SuperDARN signals which might have reached the RRI instrument.

This is a work-in-progress.


**Note**

This tool has *not* been made to be specifically portable yet. Installation 
is non-existent: you'll need to download and install DaVitPy, Numpy,
Matplotlib, and more. In addition, some of the data is retrieved from a server
from which authentication is done by the user, so unless you have access to it
already, that won't work either. Finally, it uses bash commands at a couple
different points, which is further restrictive. I apologize for all of this.

**Usage**

If you happen to have a computer with all the correct software and you have
access to the data server, the script can be used by calling it at the 
command line with the desired RRI data file as an additional argument:

If pwd is script directory (or optionally use the -i option to tell Python
to enter interactive mode after running the script):
> python script.py ./data/RRI_20160401_072744_033718_lv1_v2.h5

> python -i script.py ./data/RRI_20160401_072744_033718_lv1_v2.h5
