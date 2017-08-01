"""
plotting
--------

A module in the RRI toolkit devoted to building plots to do with the RRI
analysis, modeling, and estimation in this toolkit.



"""

import logging

# To allow us to say 'import rritk.utils', must put ../../ into path
import os, sys
rritk_dir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rritk_dir) 

import rritk.utils.data_utils as data_utils


try:
    import plots
except Exception as e:
    logging.exception("Problem importing plots: {0}".format(e))
