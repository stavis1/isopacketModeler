#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 14:51:10 2024

@author: 4vt
"""

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-m','--mzml', action = 'store', required = True,
                    help = 'The directory of .mzML files to parse.')
parser.add_argument('-p','--psms', action = 'store', required = True,
                    help = 'The proteome discoverer PSMs .txt export.')
parser.add_argument('-d','--design', action = 'store', required = True,
                    help = 'The experimental design file.')
parser.add_argument('-o','--outdir', action = 'store', required = True,
                    help = 'Directory for output files.')
parser.add_argument('-c','--cores', action = 'store', type = int, required = False, default = 1,
                    help = '[optional] Number of cores to use for mzML processing, default is 1.')
args = parser.parse_args()

from multiprocessing import Pool
from multiprocessing import Manager
import traceback
import os

import dill
import pandas as pd
import pymzml
from sortedcontainers import SortedList

from tools.fitting_tools import psm

# parse proteome discoverer PSM file
file = args.mzml[:-5]
psm_data = pd.read_csv(args.psms, sep = '\t')
psm_data = psm_data[[f.startswith(file) for f in psm_data['Spectrum File']]]
psm_cols = ['Annotated Sequence', 
            'Modifications',
            'Spectrum File', 
            'First Scan',
            'Charge', 
            'Protein Accessions']
psm_data = zip(*[psm_data[c] for c in psm_cols])

# instantiate PSM objects
psms = [psm(*d, args.label) for d in psm_data]

# parse mzML file
mzml = pymzml.run.Reader(args.mzml)
ms1s = SortedList([(s.ID,s) for s in mzml])

def process_psm(i):
    if event.is_set():
        return
    try:
        psm = psms[i]
        scan_idx = ms1s.bisect_left((psm.scan,))
        scans = ms1s[scan_idx - 3: scan_idx + 4]
        psm.parse_scans(scans)
        return psm
    except:
        traceback.print_exc()
        event.set()

#the manager contex, init_worker, and event allow the process pool to
#terminate immediately when a worker encouters a fatal error
with Manager() as manager:
    shared_event = manager.Event()
    
    def init_worker(shared_event):
        global event
        event = shared_event

    with Pool(args.cores, initializer=init_worker, initargs=(shared_event,)) as p:
        psms = p.map(process_psm, range(len(psms)))

#filter bad psms
psms = [psm for psm in psms if psm.is_useable()]

# save data
prefix = os.path.join(args.outdir, os.path.basename(args.mzml)[:-5])
with open(f'{prefix}.psms.dill', 'wb') as dillfile:
    dill.dump(psms, dillfile)

