#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 15:33:40 2024

@author: 4vt
"""

import os
from functools import cache
from copy import copy
import re
from collections import defaultdict

import numpy as np
from brainpy import isotopic_variants
from scipy.stats import binom, betabinom
from sortedcontainers import SortedList


path = os.path.abspath(__file__)[:-16]
with open(path + 'AA_masses.txt', 'r') as tsv:
    cols = tsv.readline().strip().split('\t')
    aa_formula = {l.split('\t')[0]:{e:int(c) for e,c in zip(cols[1:],l.strip().split('\t')[1:])} for l in tsv}

Δm = {'H':1.00627699,
      'C':1.00334999,
      'N':0.99704}

class hashabledict(dict):   
    def __hash__(self):
        return hash(tuple(sorted(self.items())))

    def omit(self, elm):
        newdict = self.copy()
        del newdict[elm]
        return hashabledict(newdict)

@cache
def isotope_packet(formula, charge):
    return np.array([p.intensity for p in isotopic_variants(formula, npeaks = 6, charge = charge)])

def theo_packet(psm, enrichment):
    label_dist = binom.pmf(k = range(psm['formula'][psm['label']]),
                           n = psm['formula'][psm['label']],
                           p = enrichment)
    return np.convolve(label_dist, psm['background'])

def beta_theo_packet(psm, a, b):
    label_dist = betabinom.pmf(k = range(psm['formula'][psm['label']]),
                               n = psm['formula'][psm['label']],
                               a = a,
                               b = b)
    return np.convolve(label_dist, psm['background'])

class psm:
    def __init__(self, sequence, mods, file, scan, charge, proteins, label, metadata):
        self.sequence = self.clean_seq(sequence, mods)
        self.mods = mods
        self.file = file
        self.base_name = file[:-4]
        self.scan = scan
        self.charge = charge
        self.proteins = proteins
        self.label = label
        self.metadata = metadata
        self.formula = self.calc_formula()
        self.mz = self.calc_mz()
        subformula = self.formula.omit(self.label)
        self.background = isotope_packet(subformula, self.charge)
        unenriched = isotope_packet(self.formula, self.charge)
        self.unenriched = np.concatenate((unenriched, np.zeros(len(self.mz)-len(unenriched))))
        return
    
    def clean_seq(self, seq, mods):
        seq = re.search(r'\.([^\.]+)\.',seq).group(1)
        if type(mods) == str:
            if 'N-Term(Prot)(Met-loss)' in mods:
                seq = seq[1:]
            if 'Acetyl' in mods:
                seq = '-' + seq[0].upper() + seq[1:]
        return seq
    
    def calc_mz(self):
        mz_0 = isotopic_variants(self.formula, npeaks=1, charge = self.charge)[0].mz
        init_mz = [p.mz for p in isotopic_variants(self.formula, npeaks=6, charge = self.charge)]
        mz = mz_0 + ((np.asarray(range(len(init_mz), self.formula[self.label] + 1))*Δm[self.label])/self.charge)
        comp = copy(self.formula)
        if self.label == 'C':
            comp['C[13]'] = comp.pop('C')
        elif self.label == 'N':
            comp['N[15]'] = comp.pop('N')
        terminal_mz = [p.mz for p in isotopic_variants(comp, npeaks = 30, charge = self.charge)][1:]
        mz = np.concatenate((init_mz, mz, terminal_mz), axis = None)
        return mz

    def calc_formula(self):
        formula = defaultdict(lambda: 0)
        for aa in self.sequence:
            for k,v in aa_formula[aa].items():
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
            peak_list = SortedList(zip(scan.mz, scan.i))
            peaks = [self.mymin(peak_list, mz) for mz in self.mz]
            mz_errs = np.asarray([mz - p[0] for p,mz in zip(peaks,self.mz)])
            all_errs.append(mz_errs)
            peaks = np.asarray([p[1] for p in peaks])
            all_peaks.append(peaks)
        self.intensity = np.nanmean(all_peaks, axis = 0).tolist()
        self.mz_err = np.nanmean(all_errs, axis = 0).tolist()
    
    def is_useable(self):
        return np.sum(np.isfinite(self.intensity)) > 5

class peptide:
    def __init__(self, psms):
        psm = psms[0]
        self.psms = psms
        self.sequence = psm.sequence
        self.label = psm.label
        self.formula = psm.formula
        self.mz = psm.mz
        self.background = psm.background
        self.metadata = psm.metadata
        normed = [np.asarray(p.intensity)/np.nansum(p.intensity) for p in self.psms]
        self.npeaks = max([len(n) for n in normed])
        self.obs = np.array([self.clean(np.concatenate((n, np.full(self.npeaks - len(n), np.nan)))) for n in normed])
        self.unenriched = self.reshape(psm.unenriched)
        self.mz_err = np.array([self.reshape(p.mz_err) for p in self.psms])

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
        interp = [v for pair in zip(even_interp, odd_interp) for v in pair]
        if resize:
            interp = interp[:-1]
        Δinterp = vals - interp
        Δinterp[nans] = np.full(len(vals), np.nan)[nans]
        cutoff = 0.05
        badpts = Δinterp > cutoff
        badpts[0] = False
        vals[badpts] = np.full(len(vals),np.nan)[badpts]
        return vals
    
    def reshape(self, arr):
        if len(arr) < self.npeaks:
            return np.concatenate((arr, np.zeros(self.npeaks - len(arr))))
        elif len(arr) > self.npeaks:
            return arr[:self.npeaks]
        else:
            return arr
    
    def expected(self, x):
        label = betabinom.pmf(k = range(self.formula[self.label]),
                              n = self.formula[self.label],
                              a = x[1],
                              b = x[2])
        exp = np.convolve(label, self.background)
        exp = self.reshape(exp)
        exp = (self.unenriched*(1-x[0])) + (exp*x[0])
        exp = exp/np.nansum(exp)
        return exp
    
    def interp(self, vec):
        x = np.linspace(min(self.mz),max(self.mz),256)
        vec = np.interp(x, self.mz, vec)
        return vec

    def preprocess(self):
        mz = self.mz[0]
        intensity = np.nanmean(self.obs, axis = 0)
        intensity[np.isnan(intensity)] = np.zeros(intensity.shape)[np.isnan(intensity)]
        intensity = self.interp(intensity)

        mz_std = np.nanstd(self.mz_err, axis = 0)
        mz_std = mz_std/mz
        mz_std = self.interp(mz_std)

        mz_mean = np.nanmean(self.mz_err, axis = 0)
        mz_mean = mz_mean/mz
        mz_mean = self.interp(mz_mean)
        return np.array([intensity, mz_std, mz_mean])

    def preprocess_init(self):
        intensity = np.nanmean(self.obs, axis = 0)
        intensity[np.isnan(intensity)] = np.zeros(intensity.shape)[np.isnan(intensity)]
        intensity = self.interp(intensity)
        return intensity
    




