#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 14:51:10 2024

@author: 4vt
"""

from multiprocessing import Pool, Manager
import traceback
from copy import copy
import sys

import pandas as pd
import numpy as np
import pyopenms as oms
from sortedcontainers import SortedList

from isopacketModeler.data_objects import psm, base_name, Scan

# parse PSM files into a list of data tuples
def parse_PSMs(args):
    column_names = ['sequence',
                    'file', 
                    'scan',
                    'charge', 
                    'proteinIds']
    column_map = {n:c for n,c in zip(column_names, args.PSM_headers[:len(column_names)], strict = True)}
    
    psm_data = []
    for file in args.psms:
        psm_data.append(pd.read_csv(file, sep = '\t'))
    psm_data = pd.concat(psm_data)

    #remove PSMs without matching mzML files
    bad_psms = psm_data[[base_name(f) not in args.base_names for f in psm_data[column_map['file']]]]
    psm_data = psm_data[[base_name(f) in args.base_names for f in psm_data[column_map['file']]]]
    if bad_psms.shape[0]:
        args.logs.warn('There were PSMs without a corresponding spectrum file in the design document. These will be ignored.')
        bad_files = [str(f) for f in set(bad_psms[column_map['file']])]
        args.logs.debug('The filenames for these ignored PSMs are:\n' + '\n'.join(bad_files))
    
    #add arbitrary columns listed in the optios file as a metadata dictionary
    if len(args.PSM_headers) > 5:
        psm_metadata = [{k:v for k,v in zip(d[1].keys(), d[1].values, strict = True)} for d in psm_data[args.PSM_headers[5:]].iterrows()]
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
    psm_data['design_metadata'] = [{k:v for k,v in zip(d[1].keys(), d[1].values, strict = True)} for d in metadata.iterrows()]
    psm_data['label'] = [m['label'] for m in psm_data['design_metadata']]
    psm_data['is_labeled'] = [bool(l) for l in psm_data['label']]
    
    #make duplicate control PSMs for each label used. This is for training the classifier model
    labels = sorted(set([l for l in args.design['label'] if l]))
    if not labels:
        args.logs.warning('No label elements were specified. Peptides will have both C[13] and N[15] patterns extracted. The classifier will fail.')
        labels = ['C[13]', 'N[15]']
    controls = psm_data[np.logical_not(psm_data['is_labeled'])]
    labeled = psm_data[psm_data['is_labeled']]
    psm_data = [labeled]
    for label in labels:
        temp = controls.copy()
        temp['label'] = [label]*temp.shape[0]
        psm_data.append(temp)
    psm_data = pd.concat(psm_data)
    
    # instantiate PSM objects
    psms = [psm(**d[1], args = args) for d in psm_data.iterrows()]
    args.logs.info(f'{len(psms)} PSM objects have been initialized.')
    return psms

def read_mzml(file):
    od_exp = oms.OnDiscMSExperiment()
    od_exp.openFile(file)
    ms1s = []
    for i in range(od_exp.getNrSpectra()):
        spectrum = od_exp.getSpectrum(i)
        if spectrum.getMSLevel() == 1:
            ms1s.append(Scan(spectrum))
    ms1s = SortedList(ms1s, key = lambda s: s.scan)
    return ms1s

# parse mzML files
def process_psm(psm):
    scan_idx = ms1s.bisect_left(psm)
    scans = ms1s[scan_idx - 3: scan_idx + 4]
    psm.parse_scans(scans)
    return psm if psm.is_useable() else None

def process_spectrum_data(args, psms):
    PSM_list = psms
    result_psms = []
    global ms1s
    for mzml in args.mzml_files:
        ms1s = read_mzml(mzml)
        no_extension = base_name(mzml)
        subset_psms = [p for p in PSM_list if p.base_name == no_extension]
        args.logs.debug(f'There are {len(subset_psms)} PSMs in file {no_extension}')

        with Pool(args.cores) as p:
            result_psms.extend(p for p in p.map(process_psm, subset_psms) if p is not None)
            
    args.logs.debug('Intensity data for PSMs have been extracted from mzML files.')
    args.logs.info(f'{len(result_psms)} PSMs have passed the initial usability filter.')
    return result_psms

