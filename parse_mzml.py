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
parser.add_argument('-p','--psms', action = 'append', required = True,
                    help = 'Use once per proteome discoverer PSMs .txt export.')
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

# parse proteome discoverer PSM files
psm_data = []
for file in args.psms:
    psm_data.append(pd.read_csv(file, sep = '\t'))
psm_data = pd.concat(psm_data)

mzml_files = [f for f in os.listdir(args.mzml) if f.lower().endswith('.mzml')]
good_names = set(f[:-5] for f in mzml_files)
bad_names = set(f[:-4] for f in psm_data['Spectrum File'])
if bad_names:
    print('Warning: there are files listed in the PSM data that do not have a corresponding .mzML')
    print('\n'.join(bad_names))
psm_data = psm_data[[f[:-4] in good_names for f in psm_data['Spectrum File']]]
psm_cols = ['Annotated Sequence', 
            'Modifications',
            'Spectrum File', 
            'First Scan',
            'Charge', 
            'Protein Accessions']
psm_data = zip(*[psm_data[c] for c in psm_cols])

#map filenames to metadata
design = pd.read_csv(args.design, sep = '\t')
meta_cols = [c for c in design.columns if c not in ('file', 'label')]
rows = zip(design['file'], design['label'], zip(*[design[c] for c in meta_cols]))
meta_map = {f:(l,{c:m for c,m in zip(meta_cols, meta)}) for f,l,meta in rows}

# instantiate PSM objects
psms = [psm(*d, *meta_map[d[2][:-4]]) for d in psm_data]

# parse mzML files
def process_psm(file):
    if event.is_set():
        return
    try:
        mzml = pymzml.run.Reader(f'{args.mzml}{file}')
        ms1s = SortedList([(s.ID,s) for s in mzml])
        results = []
        no_ext = file[:-5]
        subset_psms = [p for p in psms if p.file[:-4] == no_ext]
        for p in subset_psms:
            scan_idx = ms1s.bisect_left((p.scan,))
            scans = ms1s[scan_idx - 3: scan_idx + 4]
            p.parse_scans(scans)
            results.append(p)
        return results
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
        result_psms = p.map(process_psm, mzml_files)

#flatten results and filter bad psms
psms = [p for file in result_psms for p in file if p.is_useable()]

# save data
with open(f'{args.outdir}psms.dill', 'wb') as dillfile:
    dill.dump(psms, dillfile)

