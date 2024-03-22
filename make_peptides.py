#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 17:02:02 2024

@author: 4vt
"""

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-p', '--psms', action = 'store', required = True,
                    help = 'The .dill file created by parse_mzml.py')
parser.add_argument('-d', '--design', action = 'store', required = True,
                    help = 'The experimental design file')
parser.add_argument('-o','--out', action = 'store', required = True,
                    help = 'File name for results.')
args = parser.parse_args()

from collections import defaultdict
import dill
from tools.fitting_tools import peptide

with open(args.psms, 'rb') as dillfile:
    psms = dill.load(dillfile)

def fingerprint(psm):
    return (psm.sequence, psm.mods, psm.label, *psm.metadata.values())

psm_lists = defaultdict(lambda:[])
for psm in psms:
    psm_lists[fingerprint(psm)].append(psm)

peptides = [peptide(psm_list) for psm_list in psm_lists.values()]

with open(args.out, 'wb') as dillfile:
    dill.dump(peptides, dillfile)
