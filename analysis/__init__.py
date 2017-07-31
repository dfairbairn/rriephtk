"""
analysis
--------

Module for analyzing RRI experiments in terms of their transmitter-receiver
geometry or by modeling ionospheric conditions.


"""
import logging

try:
    import analysis_tools
except Exception as e:
    logging.exception("Problem importing analysis_tools: {0}".format(e))

