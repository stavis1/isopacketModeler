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
from isopacketModeler.data_objects import peptide


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

class DataGeneratingProcessTestSuite(InitializedPSMsTestSuite):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.reasonable_params = np.array([])
    
    def setUp(self):
        super().setUp()
        self.rng = np.random.default_rng(1)
        self.DGP = None
    
    def test_fit_is_poor_on_noise(self):
        for psm in self.psms:
            psm.intensity = self.rng.uniform(1e5,1e7,len(psm.mz))
            ppm5 = (psm.mz[0]/1e6)*5
            psm.mz_err = self.rng.uniform(-ppm5,ppm5,len(psm.mz))
            nan_mask = self.rng.choice(range(len(psm.mz)), 
                                       self.rng.integers(0,len(psm.mz)//3))
            psm.intensity[nan_mask] = np.nan
            psm.mz_err[nan_mask] = np.nan
        pep = peptide(self.psms)
        pep.fit_results.append(self.DGP.fit(pep))
        with self.subTest('ensure that the right number of results were added'):
            self.assertEqual(len(pep.fit_results), 1)
        with self.subTest('test that the fit is poor'):
            self.assertGreater(pep.fit_results[0].fit, 0.01)
    
    def test_fit_is_good_on_DGPs_own_expected_data(self):
        for psm in self.psms:
            psm.intensity = self.rng.uniform(1e5,1e7,len(psm.mz))
            ppm5 = (psm.mz[0]/1e6)*5
            psm.mz_err = self.rng.uniform(-ppm5,ppm5,len(psm.mz))
        pep = peptide(self.psms)
        expected_intensity = self.DGP.expected(pep, self.reasonable_params)
        expected_intensity = np.array([expected_intensity]*len(self.psms))
        pep.obs = expected_intensity
        pep.fit_results.append(self.DGP.fit(pep))
        with self.subTest('test that fit is good'):
            self.assertLess(pep.fit_results[0].fit, 0.01)
        for truth, fitted in zip(self.reasonable_params, pep.fit_results[0].params, strict = True):
            with self.subTest('test that parameters are recovered'):
                self.assertAlmostEqual(truth, fitted, delta = 0.001)
        

