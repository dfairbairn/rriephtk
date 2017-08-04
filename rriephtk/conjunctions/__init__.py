"""
conjunct
--------

Module for analyzing conjunction of ePOP-RRI experiments with SuperDARN 
transmissions and data.


"""
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

# To allow us to say 'import rritk.utils', must put ../../ into path
import os, sys
rritk_dir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rritk_dir) 

import rritk.utils.data_utils as data_utils

try:
    import conjunctions
except Exception as e:
    logging.exception("Problem importing time_align: {0}".format(e))

try:
    import time_align
except Exception as e:
    logging.exception("Problem importing time_align: {0}".format(e))
