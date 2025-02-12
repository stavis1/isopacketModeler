#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 17:02:02 2024

@author: 4vt
"""

from collections import defaultdict
from isopacketModeler.data_objects import peptide

def fingerprint(psm):
    return (psm.raw_sequence, psm.file)

def initialize_peptides(args, psms, bad_psms):
    psm_lists = defaultdict(lambda:[])
    for psm in psms:
        psm_lists[fingerprint(psm)].append(psm)
    
    good_fingerprints = set(psm_lists.keys())
    for psm in bad_psms:
        this_fingerprint = fingerprint(psm)
        if this_fingerprint in good_fingerprints:
            psm_lists[this_fingerprint].append(psm)
    
    peptides = [peptide(psm_list) for psm_list in psm_lists.values()]
    args.logs.info(f'{len(peptides)} Peptides have been identified.')
    return peptides
