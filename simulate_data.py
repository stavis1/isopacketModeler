#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 16:42:24 2024

@author: 4vt
"""

from tools.options import options
args = options.alt_init('--options options.toml')

import numpy as np
from scipy.stats import binom
from brainpy import isotopic_variants
import dill

rng = np.random.default_rng(1)

#collect control peptides
fingerprint_cols =  [c for c in args.design.columns if not c in ['file','label', 'control']]
fingerprints = zip(*[args.design[c] for c in fingerprint_cols])
controls = set(f for f,c in zip(fingerprints,args.design['control']) if c)

def fingerprint(pep):
    return tuple(pep.metadata[c] for c in fingerprint_cols)

with open('peptides.dill','rb') as dillfile:
    peptides = dill.load(dillfile)
ctrl_peps = [p for p in peptides if fingerprint(p) in controls]

#set up control data
features = np.array([p.preprocess_init() for p in ctrl_peps])
labels = np.zeros(features.shape[0])

#initialize noise

#simulate expected data

#format outputs



