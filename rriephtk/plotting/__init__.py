"""
plotting
--------

A module in the RRI toolkit devoted to building plots to do with the RRI
analysis, modeling, and estimation in this toolkit.



"""

import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

# To allow us to say 'import rriephtk.utils', must put ../../ into path
import os, sys
rriephtk_dir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rriephtk_dir) 

import rriephtk.utils.data_utils as data_utils


try:
    import plots
except Exception as e:
    logging.exception("Problem importing plots: {0}".format(e))
