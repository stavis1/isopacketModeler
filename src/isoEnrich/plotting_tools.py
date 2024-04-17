#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 16:50:17 2024

@author: 4vt
"""


import matplotlib.pyplot as plt
from matplotlib import cm

import numpy as np

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})


def get_colors(vals):
    low = min(vals)
    high = max(vals)
    return [cm.plasma(int(((val-low)/(high-low))*cm.plasma.N)) for val in vals]

def get_sm(vals):
    colormap = plt.cm.get_cmap('plasma')
    sm = plt.cm.ScalarMappable(cmap=colormap)
    sm.set_clim(vmin = min(vals), vmax = max(vals))
    return sm
