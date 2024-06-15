#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 17:02:02 2024

@author: 4vt
"""

from isoEnrich.options import options
args = options.alt_init()

from collections import defaultdict
from tools.fitting_tools import peptide

def fingerprint(psm):
    return (psm.sequence, psm.mods, psm.label, *psm.metadata.values())

def initialize_peptides(psms):
    psm_lists = defaultdict(lambda:[])
    for psm in psms:
        psm_lists[fingerprint(psm)].append(psm)
    
    peptides = [peptide(psm_list) for psm_list in psm_lists.values()]
    return peptides
