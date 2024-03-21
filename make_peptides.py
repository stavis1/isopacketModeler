#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 17:02:02 2024

@author: 4vt
"""

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-p', '--psms', action = 'append', required = True,
                    help = 'Use once per .dill file created by parse_mzml.py')
parser.add_argument('-o','--out', action = 'store', required = True,
                    help = 'File name for results.')
args = parser.parse_args()

from collections import defaultdict
import dill
from tools.fitting_tools import peptide

labels = []
psms = []
for file in args.psms:
    with open(file, 'rb') as dillfile:
        psms.extend(dill.load(dillfile))
    labels.append(psms[-1].label)

labels = set(labels)
if len(labels) > 1:
    raise Exception('Multiple labels found. This script should only be used to combine replicates.')

def fingerprint(psm):
    return (psm.sequence, psm.mods)

psm_lists = defaultdict(lambda:[])
for psm in psms:
    psm_lists[fingerprint(psm)].append(psm)

peptides = [peptide(psm_list) for psm_list in psm_lists]

with open(args.out, 'wb') as dillfile:
    dill.dump(peptides, dillfile)
