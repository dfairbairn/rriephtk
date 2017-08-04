"""
analysis
--------

Module for analyzing RRI experiments in terms of their transmitter-receiver
geometry or by modeling ionospheric conditions.


"""
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

# To allow us to say 'import rriephtk.utils', must put ../../ into path
import os, sys
rriephtk_dir = os.path.dirname(os.path.dirname(os.getcwd()))
#sys.path.insert(0, rritk_dir) 
sys.path.append(rriephtk_dir)

import rriephtk.utils.data_utils as data_utils


try:
    import analysis_tools
except Exception as e:
    logging.exception("Problem importing analysis_tools: {0}".format(e))

try:
    import magnet_data
except Exception as e:
    logging.exception("Problem impoting magnet_data: {0}".format(e))
