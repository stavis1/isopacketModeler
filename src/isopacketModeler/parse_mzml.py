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

from isopacketModeler.data_objects import psm, base_name

# parse PSM files into a list of data tuples
def parse_PSMs(args):
    column_names = ['sequence',
                    'file', 
                    'scan',
                    'charge', 
                    'proteinIds']
    column_map = {n:c for n,c in zip(column_names, args.PSM_headers[:len(column_names)])}
    
    psm_data = []
    for file in args.psms:
        psm_data.append(pd.read_csv(file, sep = '\t'))
    psm_data = pd.concat(psm_data)

    #remove PSMs without matching mzML files
    ninit = psm_data.shape[0]
    psm_data = psm_data[[base_name(f) in args.base_names for f in psm_data[column_map['file']]]]
    if psm_data.shape[0] < ninit:
        args.logs.warn('There were PSMs without a corresponding spectrum file in the design document. These will be ignored.')
    
    #add arbitrary columns listed in the optios file as a metadata dictionary
    if len(args.PSM_headers) > 5:
        psm_metadata = [{k:v for k,v in zip(d[1].keys(), d[1].values)} for d in psm_data[args.PSM_headers[5:]].iterrows()]
    else:
        psm_metadata = [{}]*psm_data.shape[0]
    psm_data['psm_metadata'] = psm_metadata
    
    #subset the dataframe to only the columns we use
    psm_data = psm_data[args.PSM_headers[:5] + ['psm_metadata']]
    psm_data.columns = ['raw_sequence', 'file', 'scan', 'charge', 'proteins', 'psm_metadata']
    return psm_data

def initialize_psms(args, psm_data):
    #gather initialization data for PSMs
    design_data = copy(args.design)
    design_data.index = [base_name(f) for f in design_data['file']]

    #add design metadata dictionaries to PSMs    
    metadata = design_data.loc[[base_name(f) for f in psm_data['file']]]
    psm_data['design_metadata'] = [{k:v for k,v in zip(d[1].keys(), d[1].values)} for d in metadata.iterrows()]
    psm_data['label'] = [m['label'] for m in psm_data['design_metadata']]
    
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
    
    # instantiate PSM objects
    psms = [psm(**d[1], args = args) for d in psm_data.iterrows()]
    args.logs.info('PSM objects have been initialized.')
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
    logger.info(f'Raw isotope patterns have been extracted from {file}.')
    return results

def process_spectrum_data(args, psms):
    global PSM_list
    PSM_list = psms
    global logger
    logger = args.logs

    with Pool(args.cores) as p:
        result_psms = p.map(process_psm, args.mzml_files)

    #flatten results and filter bad psms
    psms = [p for file in result_psms for p in file if p.is_useable()]
    return psms

