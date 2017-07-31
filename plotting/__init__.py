"""
plotting
--------

A module in the RRI toolkit devoted to building plots to do with the RRI
analysis, modeling, and estimation in this toolkit.



"""

import logging

try:
    import fig1_plots
except Exception as e:
    logging.exception("Problem importing fig1_plots: {0}".format(e))
