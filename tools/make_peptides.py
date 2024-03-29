#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 17:02:02 2024

@author: 4vt
"""

from tools.options import options
args = options.alt_init()

from collections import defaultdict
import dill
from tools.fitting_tools import peptide

with open('psms.dill', 'rb') as dillfile:
    psms = dill.load(dillfile)

def fingerprint(psm):
    return (psm.sequence, psm.mods, psm.label, *psm.metadata.values())

psm_lists = defaultdict(lambda:[])
for psm in psms:
    psm_lists[fingerprint(psm)].append(psm)

peptides = [peptide(psm_list) for psm_list in psm_lists.values()]

with open('peptides.dill', 'wb') as dillfile:
    dill.dump(peptides, dillfile)
