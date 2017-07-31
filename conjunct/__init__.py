"""
conjunct
--------

Module for analyzing conjunction of ePOP-RRI experiments with SuperDARN 
transmissions and data.


"""
import logging

try:
    import conjunctions
except Exception as e:
    logging.exception("Problem importing time_align: {0}".format(e))

try:
    import time_align
except Exception as e:
    logging.exception("Problem importing time_align: {0}".format(e))
