"""
file: 'conjunctions.py'
description:
    This is the 'script 2.0' file. The functions herein allow 
    determination of intersections in CASSIOPE ephemeris during which
    it passes through the field-of-view of a SuperDARN radar.

    Basic functions simply report the intersection, and additionally
    the .h5 data files can have ephemeris points tagged with the 
    intersection details, and sections from the SuperDARN data logs can
    even be retrieved.

author: David Fairbairn
date: May 2017

This is a reworking/improvement of the original 'script.py' from 2016.
"""

import os
import subprocess
import sys
import logging

import davitpy
from davitpy import pydarn
import timeit
import math

import datetime as dt
import numpy as np

import data_utils


