"""
This file simply can be imported by example scripts in this folder in order to
add the rriephtk root directory to the python path.
"""

from inspect import getsourcefile
import os.path
import sys

cur_path = os.path.abspath(getsourcefile(lambda:0))
cur_dir = os.path.dirname(cur_path)
parent_dir = cur_dir[:cur_dir.rfind(os.path.sep)]

sys.path.insert(0, parent_dir)
