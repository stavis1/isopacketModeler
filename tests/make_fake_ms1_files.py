#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 14:41:03 2024

@author: 4vt
"""

from brainpy import isotopic_variants
import numpy as np
import pandas as pd
rng = np.random.default_rng()

header = '''H	CreationDate Thu Apr 18 19:36:12 2024
H	Extractor	ProteoWizard
H	Extractor version	Xcalibur
H	Source file	SLT10_QEDIA_firstmassTest.raw
'''


class scan:
    def __init__(self, ms_level, enrichment, scan, formula = None):
        self.ms_level = ms_level
        self.enrichment = enrichment
        self.scan = scan
        self.rt = scan/10
        if ms_level == 2:
            self.parse_formula(formula)
        if ms_level == 2 or enrichment < 0:
            self.generate_random()
        else:
            self.generate_enriched(self.enrichment)

    def parse_formula(self, formula):
        isotopes = isotopic_variants(formula, npeaks = 20, charge = 2)
        self.precursor = isotopes[0].intensity
        isotope_Δ_mz = 1.00334999/2
        mz = [i.mz for i in isotopes]
        # mz += 
        

    def generate_random(self):
        mz = sorted(rng.uniform(400, 2000, 50))
        intensity = rng.uniform(1e5, 1e7, 50)
        self.spectrum = list(zip(mz, intensity))
    
    def __str__(self):
        base_peak = max(self.spectrum, key = lambda x: x[1])
        #write scan number and optionally precursor m/z
        if self.ms_level == 2:
            entry = f'S\t{self.scan}\t{self.scan}\t{self.precursor}\n'
        else:
            entry = f'S\t{self.scan}\t{self.scan}\n'
        #write retention time
        entry += f'I\tNativeID\tcontrollerType=0 controllerNumber=1 scan={self.scan}\n'
        entry += f'I\tRTime\t{self.rt}\n'
        entry += f'I\tBPI\t{base_peak[1]}\n'
        entry += f'I\tBPM\t{base_peak[0]}\n'
        if self.ms_level == 2:
            entry += f'Z\t2\t{self.precursor*2}\n'
        #write spectrum peaks
        entry += '\n'.join([f'{mz} {i}' for mz, i in self.spectrum]) + '\n'
        return entry


scans = []
for i in range(1, 20):
    scans.append(scan(1,-1,i))

def write_files(prefix, scans):
    with open(prefix + '.ms1', 'w') as ms1:
        ms1.write(header)
        for scan in scans:
            if scan.ms_level == 1:
                ms1.write(str(scan))

write_files('test', scans)

