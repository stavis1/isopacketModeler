#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 14:51:10 2024

@author: 4vt
"""

from isoEnrich.options import options
args = options()

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
def parse_PD(args):
    psm_data = []
    for file in args.psms:
        psm_data.append(pd.read_csv(file, sep = '\t'))
    psm_data = pd.concat(psm_data)
    psm_data = psm_data[[f[:-4] in args.base_names for f in psm_data['Spectrum File']]]
    psm_cols = ['Annotated Sequence', 
                'Modifications',
                'Spectrum File', 
                'First Scan',
                'Charge', 
                'Protein Accessions']
    psm_data = zip(*[psm_data[c] for c in psm_cols])
    return psm_data

def initialize_psms(args, psm_data):
    #gather initialization data for PSMs
    meta_cols = [c for c in args.design.columns if c not in ('file', 'label')]
    meta_rows = zip(*[args.design[c] for c in meta_cols])
    rows = list(zip(args.design['file'], args.design['label'], meta_rows))
    
    #make duplicate control PSMs for each label used. This is for training the classifier model
    ctrl_idxs = [i for i,r in enumerate(rows) if not r[1]]
    ctrls = [rows[i] for i in ctrl_idxs]
    rows = [r for i,r in enumerate(rows) if not i in ctrl_idxs]
    labels = set(l for l in args.design['label'] if l)
    rows.extend([(r[0],l,r[2]) for l in labels for r in ctrls])
    
    meta_map = {f:(l,{c:m for c,m in zip(meta_cols, meta)}) for f,l,meta in rows}
    
    # instantiate PSM objects
    psms = [psm(*d, *meta_map[d[2][:-4]]) for d in psm_data]
    return psms

# parse mzML files
def process_psm(file):
    if event.is_set():
        return
    try:
        mzml = pymzml.run.Reader(file)
        ms1s = SortedList([(s.ID,s) for s in mzml])
        results = []
        no_ext = os.path.basename(file)[:-5]
        subset_psms = [p for p in PSM_list if p.base_name == no_ext]
        for p in subset_psms:
            scan_idx = ms1s.bisect_left((p.scan,))
            scans = ms1s[scan_idx - 3: scan_idx + 4]
            p.parse_scans(scans)
            results.append(p)
        return results
    except:
        traceback.print_exc()
        event.set()

def process_spectrum_data(args, psms):
    global PSM_list
    PSM_list = psms
    #the manager contex, init_worker, and event allow the process pool to
    #terminate immediately when a worker encouters a fatal error
    with Manager() as manager:
        shared_event = manager.Event()
        
        def init_worker(shared_event):
            global event
            event = shared_event
            
        
        with Pool(args.cores, initializer=init_worker, initargs=(shared_event,)) as p:
            result_psms = p.map(process_psm, args.mzml_files)

    #flatten results and filter bad psms
    psms = [p for file in result_psms for p in file if p.is_useable()]
    return psms

# def save_psms(psms):
#     # save data
#     with open('psms.dill', 'wb') as dillfile:
#         dill.dump(psms, dillfile)

