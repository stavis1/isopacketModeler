#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 15:33:40 2024

@author: 4vt
"""

from functools import cache
from copy import copy
import re
import os
from collections import defaultdict

import numpy as np
from brainpy import isotopic_variants
from sortedcontainers import SortedList


Δm = {'H[2]':1.0062767458900002,
      'C[13]':1.0033548350700006,
      'N[15]':0.9970348944500014,
      'O[17]':1.0042171369299986,
      'O[18]':2.0042449932900013,
      'S[33]':0.9993877353999991,
      'S[34]':1.9957958295999987,
      'S[36]':3.9950095355999977}

class hashabledict(dict):   
    def __hash__(self):
        return hash(tuple(sorted(self.items())))

    def omit(self, elm):
        newdict = self.copy()
        del newdict[elm]
        return hashabledict(newdict)

def base_name(filename):
    return os.path.basename(os.path.splitext(filename)[0])

@cache
def isotope_packet(formula, charge):
    return np.array([p.intensity for p in isotopic_variants(formula, npeaks = 6, charge = charge)])

class psm:
    def __init__(self,
                 raw_sequence,
                 file, 
                 scan,
                 charge,
                 proteins,
                 psm_metadata,
                 design_metadata,
                 label,
                 is_labeled,
                 args):
        
        self.raw_sequence = raw_sequence
        self.file = file
        self.scan = scan
        self.charge = charge
        self.proteins = proteins
        self.psm_metadata = psm_metadata
        self.design_metadata = design_metadata
        self.label = label
        self.label_elm = re.search(r'[A-Z][a-z]?', label).group()
        self.is_labeled = is_labeled
        self.AA_formulae = args.AA_formulae
        
        self.sequence = self.clean_seq(self.raw_sequence)
        self.base_name = base_name(self.file)

        with open(self.AA_formulae, 'r') as tsv:
            cols = tsv.readline().strip().split('\t')
            self.aa_formulae = {l.split('\t')[0]:{e:int(c) for e,c in zip(cols[1:],l.strip().split('\t')[1:], strict = True)} for l in tsv}
        self.formula = self.calc_formula()
        self.mz = self.calc_mz()
        subformula = self.formula.omit(self.label_elm)
        self.background = isotope_packet(subformula, self.charge)
        unenriched = isotope_packet(self.formula, self.charge)
        self.unenriched = np.concatenate((unenriched, np.zeros(len(self.mz)-len(unenriched))))
        return
    
    def clean_seq(self, seq):
        #The regex excludes non bracket characters at the beginning or end of the string that are demarcated from 
        #the middle with a period, e.g. T.ES.T -> ES, this allows it to work with Proteome discoverer annotated sequences.
        #The capturing group in the middle grabs a sequence of residues. These residues start with a character that is
        #neither whitespace nor an open bracket and are optionally followed by an arbitrary number of characters contained
        #by brackets, the only constraint on these characters is that they do not contain a close bracket.
        seq = re.search(r'\A(?:[^\.\[]+\.)?((?:[^\[\s\.](?:\[[^\]]+\])?)+)(?:\.[^\.\]]+)?\Z',seq).group(1)
        return seq
    
    def calc_mz(self):
        mz_0 = isotopic_variants(self.formula, npeaks=1, charge = self.charge)[0].mz
        init_mz = [p.mz for p in isotopic_variants(self.formula, npeaks=6, charge = self.charge)]
        mz = mz_0 + ((np.asarray(range(len(init_mz), self.formula[self.label_elm] + 1))*Δm[self.label])/self.charge)
        comp = copy(self.formula)
        comp[self.label] = comp.pop(self.label_elm)
        terminal_mz = [p.mz for p in isotopic_variants(comp, npeaks = 30, charge = self.charge)][1:]
        mz = np.concatenate((init_mz, mz, terminal_mz), axis = None)
        return mz

    def calc_formula(self):
        formula = defaultdict(lambda: 0)
        residues = re.findall(r'[^\[](?:\[[^\]]+\])?', self.sequence)
        for aa in residues:
            for k,v in self.aa_formulae[aa].items():
                formula[k] += v
        formula['H'] += 2
        formula['O'] += 1
        formula = hashabledict(formula)
        return formula
    
    def mymin(self, peaks, mz):
        idx = max(1,peaks.bisect_left((mz,)))
        query = min(peaks[idx-1:idx+1], 
                    key = lambda x: abs(x[0] - mz), 
                    default = (np.nan,))
        if abs(query[0] - mz)  < (10/1e6)*mz:
            return query
        else:
            return (np.nan, np.nan)
    
    def parse_scans(self, scans):
        all_peaks = []
        all_errs = []
        for i,scan in scans:
            peak_list = SortedList(zip(scan.get_mz_array(), scan.get_intensity_array(), strict = True))
            peaks = [self.mymin(peak_list, mz) for mz in self.mz]
            mz_errs = np.asarray([mz - p[0] for p,mz in zip(peaks,self.mz, strict = True)])
            all_errs.append(mz_errs)
            peaks = np.asarray([p[1] for p in peaks])
            all_peaks.append(peaks)
        self.intensity = np.nanmean(all_peaks, axis = 0).tolist()
        self.mz_err = np.nanmean(all_errs, axis = 0).tolist()
    
    def is_useable(self):
        finite = np.isfinite(self.intensity)
        if np.sum(finite) > 5:
            run = 0
            for i in finite:
                if i:
                    run += 1
                else:
                    run = 0
                if run > 4:
                    return True
        return False

class peptide:
    def __init__(self, psms):
        psm = psms[0]
        self.psms = psms
        self.sequence = psm.sequence
        self.raw_sequence = psm.raw_sequence
        self.label = psm.label
        self.label_elm = psm.label_elm
        self.formula = psm.formula
        self.mz = psm.mz
        self.background = psm.background
        self.design_metadata = psm.design_metadata
        self.psm_metadata = [p.psm_metadata for p in psms]
        normed = [np.asarray(p.intensity)/np.nansum(p.intensity) for p in self.psms]
        self.npeaks = max([len(n) for n in normed])
        self.obs = np.array([self.clean(np.concatenate((n, np.full(self.npeaks - len(n), np.nan)))) for n in normed])
        self.unenriched = self.reshape(psm.unenriched)
        self.mz_err = np.array([self.reshape(p.mz_err) for p in self.psms])
        self.fit_results = []
        self.canonical_fit = None

    def clean(self, vals):
        vals = copy(vals)
        #remove singletons
        for i in range(1, len(vals)-1):
            if np.sum(np.isnan(vals[i-1:i+2])) > 1:
                vals[i] = np.nan
        #remove points far from the trend
        nans = np.isnan(vals)
        vals[nans] = np.zeros(len(vals))[nans]
        evens = vals[np.array(range(len(vals)))%2 == 0]
        odds = vals[np.array(range(len(vals)))%2 == 1]
        odd_interp = np.interp(range(1,len(vals), 2), range(0,len(vals), 2), evens)
        even_interp = np.interp(range(0,len(vals), 2), range(1,len(vals), 2), odds)
        if len(odd_interp) != len(even_interp):
            resize = True
            odd_interp = np.concatenate((odd_interp,[np.nan]))
        else:
            resize = False
        interp = [v for pair in zip(even_interp, odd_interp, strict = True) for v in pair]
        if resize:
            interp = interp[:-1]
        Δinterp = vals - interp
        Δinterp[nans] = np.full(len(vals), np.nan)[nans]
        cutoff = 0.05
        badpts = Δinterp > cutoff
        badpts[0] = False
        vals[badpts] = np.full(len(vals),np.nan)[badpts]
        vals = vals/np.nansum(vals)
        return vals
    
    def reshape(self, arr):
        if len(arr) < self.npeaks:
            return np.concatenate((arr, np.zeros(self.npeaks - len(arr))))
        elif len(arr) > self.npeaks:
            return arr[:self.npeaks]
        else:
            return arr
    
    def report(self):
        def good_elm(elm):
            return type(elm) == 'str' or not hasattr(elm, '__iter__')
        
        fields = defaultdict(lambda:None)
        fields['peptide'] = self.psms[0].raw_sequence
        fields['proteins'] = self.psms[0].proteins
        fields.update({k:v for k,v in self.__dict__.items() if good_elm(v)})
        fields.update(self.design_metadata)
        fields['PSM_count'] = len(self.psm_metadata)
        psm_metadata_cols = set(k for m in self.psm_metadata for k in m.keys())
        for metadata_col in psm_metadata_cols:
            vals = [m[metadata_col] for m in self.psm_metadata]
            if np.all(np.isreal(vals)):
                fields[f'PSMs_mean_{metadata_col}'] = np.mean(vals)
            else:
                fields[f'PSMs_unique_{metadata_col}'] = ';'.join([str(v) for v in set(vals)])
        fields['canonical_DGP'] = self.canonical_fit.dgp_name
        fields.update({f'canonical_{k}':v for k,v in self.canonical_fit.__dict__.items() if good_elm(v)})
        fields.update({f'canonical_param{i}':v for i,v in enumerate(self.canonical_fit.params)})
        for dgp in self.fit_results:
            fields.update({f'{dgp.dgp_name}_{k}':v for k,v in dgp.__dict__.items() if good_elm(v)})
            fields.update({f'{dgp.dgp_name}_param{i}':v for i,v in enumerate(dgp.params)})           
        return fields


