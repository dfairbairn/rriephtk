#import os
#import ConfigParser as cp
#parsing config files
#config = core.parse_config_file()
#hdw_files_path = config.get("core","hdw_files_path")
#hdw_info = core.parse_hdw_files(hdw_files_path)

import logging

__all__ = ["utils", "analysis", "plotting", "conjunctions"]

try:
    from . import utils
except Exception as e:
    logging.exception("Problem importing utils: {0}".format(e))

try:
    from . import conjunctions
except Exception as e:
    logging.exception("Problem importing conjunctions: {0}".format(e))

try:
    from . import analysis
except Exception as e:
    logging.exception("Problem importing analysis: {0}".format(e))

try:
    from . import plotting 
except Exception as e:
    logging.exception("Problem importing plotting: {0}".format(e))
