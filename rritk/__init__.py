"""
rritk
-----

A toolkit for analyzing RRI experiments.

"""
import logging

try:
    import utils
except Exception as e:
    logging.exception("Problem importing utils: {0}".format(e))

try:
    import conjunctions 
except Exception as e:
    logging.exception("Problem importing conjunctions: {0}".format(e))

try:
    import analysis 
except Exception as e:
    logging.exception("Problem importing analysis: {0}".format(e))

try:
    import plotting 
except Exception as e:
    logging.exception("Problem importing plotting: {0}".format(e))
