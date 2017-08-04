"""
utils
-----

Module for miscellaneous small utilities that are used by this toolkit
but don't belong anywhere specifically


"""
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

try:
    import data_utils
except Exception as e:
    logging.exception("Problem importing data_utils: {0}".format(e))

try:
    import script_utils
except Exception as e:
    logging.exception("Problem importing data_utils: {0}".format(e))

