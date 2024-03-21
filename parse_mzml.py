#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 14:51:10 2024

@author: 4vt
"""

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-m','--mzml', action = 'store', required = True,
                    help = 'The .mzML file to parse.')
parser.add_argument('-p','--psms', action = 'store', required = True,
                    help = 'The proteome discoverer PSMs .txt export.')
parser.add_argument('-l','--label', action = 'store', required = True,
                    help = 'The element symbol used for stable isotope labeling, i.e. C for 13C labeling.')
parser.add_argument('-o','--outdir', action = 'store', required = True,
                    help = 'Directory for output files.')
parser.add_argument('-c','--cores', action = 'store', type = int, required = False, default = 1,
                    help = '[optional] Number of cores to use for mzML processing, default is 1.')
args = parser.parse_args()

from multiprocessing import Pool

import dill
import numpy as np
import pandas as pd
import pymzml
from sortedcontainers import SortedList

from tools.fitting_tools import psm

# parse proteome discoverer PSM file

# instantiate PSM objects

# parse mzML file

# save data


