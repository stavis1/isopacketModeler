#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 14:51:10 2024

@author: 4vt
"""

from multiprocessing import Pool
import os
from copy import copy

import pandas as pd
import pymzml
from sortedcontainers import SortedList

from isopacketModeler.fitting_tools import psm, base_name

# parse PSM files into a list of data tuples
def parse_PSMs(args):
    psm_data = []
    for file in args.psms:
        psm_data.append(pd.read_csv(file, sep = '\t'))
    psm_data = pd.concat(psm_data)
    ninit = psm_data.shape[0]
    psm_data = psm_data[[base_name(f) in args.base_names for f in psm_data['Spectrum File']]]
    if psm_data.shape[0] < ninit:
        args.logs.warn('There were PSMs without a corresponding spectrum file in the design document. These will be ignored.')
    psm_data = psm_data[args.PSM_headers]
    psm_data.columns = ['seq', 'file', 'scan', 'charge', 'proteins'] + args.PSM_headers[5:]
    return psm_data

def initialize_psms(args, psm_data):
    #gather initialization data for PSMs
    design_data = copy(args.design)
    design_data.index = [base_name(f) for f in design_data['file']]

    #add metadata to PSMs    
    metadata = design_data.loc[[base_name(f) for f in psm_data['file']]]
    metadata.index = range(metadata.shape[0])
    psm_data = pd.concat((psm_data, metadata), axis = 1)    
    
    #make duplicate control PSMs for each label used. This is for training the classifier model
    labels = sorted(set([l for l in args.design['label'] if l]))
    controls = psm_data[psm_data['label'] == '']
    labeled = psm_data[psm_data['label'] != '']
    psm_data = [labeled]
    for label in labels:
        temp = controls.copy()
        temp['label'] = [label]*temp.shape[0]
        psm_data.append(temp)
    psm_data = pd.concat(psm_data)
    psm_data['aa_formulae'] = [args.AA_formulae]*psm_data.shape[0]
    
    # instantiate PSM objects
    psms = [psm({k:v for k,v in zip(d[1].keys(), d[1].values)}) for d in psm_data.iterrows()]
    return psms

def read_mzml(file):
    mzml = pymzml.run.Reader(file)
    ms1s = SortedList([(s.ID,s) for s in mzml if s.ms_level == 1])
    mzml.close()
    return ms1s

# parse mzML files
def process_psm(file):
    ms1s = read_mzml(file)
    results = []
    no_extension = os.path.basename(file)[:-5]
    subset_psms = [p for p in PSM_list if p.base_name == no_extension]
    for p in subset_psms:
        scan_idx = ms1s.bisect_left((p.scan,))
        scans = ms1s[scan_idx - 3: scan_idx + 4]
        p.parse_scans(scans)
        results.append(p)
    return results

def process_spectrum_data(args, psms):
    global PSM_list
    PSM_list = psms

    with Pool(args.cores) as p:
        result_psms = p.map(process_psm, args.mzml_files)

    #flatten results and filter bad psms
    psms = [p for file in result_psms for p in file if p.is_useable()]
    return psms

