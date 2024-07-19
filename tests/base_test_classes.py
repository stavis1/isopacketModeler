#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 17:52:16 2024

@author: 4vt
"""
import sys
import os
import unittest
import shutil
import warnings

import numpy as np
import pandas as pd
from brainpy import isotopic_variants
from sortedcontainers import SortedList

from isopacketModeler.options import options
from isopacketModeler.parse_mzml import initialize_psms, process_psm, process_spectrum_data, read_mzml


class ParsedOptionsTestSuite(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.init_argv = sys.argv
        
    def setUp(self):
        sys.argv += '--options options.toml'.split()
        self.args = options()
    
    def tearDown(self):
        sys.argv = self.init_argv
        shutil.rmtree('temp')
        if os.path.exists('tests.log'):
            shutil.rmtree('tests.log')

class InitializedPSMsTestSuite(ParsedOptionsTestSuite):
    def setUp(self):
        super().setUp()
        self.N = 3
        self.raw_sequence = ['TEST', 'A.TEST.A', 'TES[1.23]T']
        self.file = ['test1.mzML']*self.N
        self.scan = range(self.N)
        self.charge = [2]*self.N
        self.prots = ['TESTPROT']*self.N
        self.psm_data = pd.DataFrame({'raw_sequence':self.raw_sequence,
                                      'file':self.file,
                                      'scan':self.scan,
                                      'charge':self.charge,
                                      'proteins':self.prots,
                                      'psm_metadata':[{}]*self.N})
        self.psms = initialize_psms(self.args, self.psm_data)
        
class ProcessedPSMsTestSuite(InitializedPSMsTestSuite):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.read_mzml = read_mzml
        self.process_psm = process_psm
    
    def setUp(self):
        super().setUp()
        #initialize PSM list    
        PSM_list = self.psms
        for psm in PSM_list:
            psm.scan = 10
        #simulate unlabeled spectrum
        spectrum = isotopic_variants(PSM_list[0].formula, npeaks = 6, charge = PSM_list[0].charge)
        mz = [p.mz for p in spectrum]
        intensity = [p.intensity for p in spectrum]
        #construct ms1 list
        class mock_ms1():
            def __init__(self):
                self.i = intensity
                self.mz = mz
        ms1s = SortedList([(i, mock_ms1()) for i in range(20)])

        def read_mzml(file):
            return ms1s
        
        #inject our mocked data into the namespace of the function
        process_psm.__globals__['read_mzml'] = read_mzml
        process_psm.__globals__['PSM_list'] = PSM_list
        process_spectrum_data.__globals__['process_psm'] = process_psm
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.psms = process_spectrum_data(self.args, PSM_list)
            assert len(self.psms) == self.N

    def tearDown(self):
        super().tearDown()
        globals()['process_psm'] = self.process_psm
        globals()['read_mzml'] = self.read_mzml

class DataGeneratingProcessTestSuite(ParsedOptionsTestSuite):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.DGP = None
    
    def setUp(self):
        super().setUp()
        self.rng = np.random.default_rng(1)
    
    def test_fit_is_poor_on_noise(self):
        pass
    
    def test_fit_is_good_on_DGPs_own_expected_data(self):
        pass
    
    
    

