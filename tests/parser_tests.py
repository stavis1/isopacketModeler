#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 18:14:31 2024

@author: anon
"""

import unittest
import warnings

from brainpy import isotopic_variants
from sortedcontainers import SortedList
import numpy as np

import base_test_classes
from isopacketModeler.parse_mzml import parse_PSMs, initialize_psms, process_psm, process_spectrum_data, read_mzml

class parsePSMsTestSuite(base_test_classes.ParsedOptionsTestSuite):
    def test_parses_PD_psms_file(self):
        psms = parse_PSMs(self.args)
        with self.subTest('Test that 6 psms were found'):
            self.assertEqual(psms.shape[0], 10)
        for psm in psms:
            with self.subTest('Test that PSM have the right number of columns'):
                self.assertEqual(psms.shape[1], 6)
    
    def test_parser_collects_arbitrary_metadata(self):
        meta_keys = ["PSMs Workflow ID", "Tags"]
        self.args.PSM_headers += meta_keys
        psms = parse_PSMs(self.args)
        for meta in psms['psm_metadata']:
            with self.subTest('Test that metadata dictionaries have the right entries'):
                self.assertSequenceEqual(list(meta.keys()), meta_keys)

class initializePSMsTestSuite(base_test_classes.InitializedPSMsTestSuite):
    def test_label_file_PSMs_initialize_correctly(self):
        N = self.N
        psms = self.psms
        with self.subTest('Check the number of PSMS'):
            self.assertEqual(len(psms), N)
        with self.subTest('Check sequence cleaning is correct'):
            self.assertListEqual([p.sequence for p in psms], ['TEST', 'TEST', 'TES[1.23]T'])
        with self.subTest('Check the scan numbers are corret'):
            self.assertListEqual([p.scan for p in psms], list(range(N)))
        with self.subTest('Check the charges are corret'):
            self.assertListEqual([p.charge for p in psms], [2]*N)
        with self.subTest('Check the proteins are corret'):
            self.assertListEqual([p.proteins for p in psms], ['TESTPROT']*N)
        with self.subTest('Check the labels are corret'):
            self.assertListEqual([p.label for p in psms], ['C[13]']*N)
        
        ref_formula = {'C':16, 'H':28, 'N':4, 'O':10, 'S':0, 'Se':0}
        for elm in ref_formula.keys():
            with self.subTest('test formula is correct'):
                self.assertEqual(ref_formula[elm], psms[0].formula[elm])
    
    def test_control_PSMs_are_duplicated(self):
        psm_data = self.psm_data
        psm_data['file'] = ['test2.mzML']*self.N
        psms = initialize_psms(self.args, psm_data)
        with self.subTest('Check the number of PSMS'):
            self.assertEqual(len(psms), self.N*2)
        with self.subTest('Check both labels are used'):
            self.assertListEqual([p.label for p in psms], ['C[13]']*self.N + ['N[15]']*self.N)

class mzmlReaderTestSuite(base_test_classes.ParsedOptionsTestSuite):
    def test_mzml_reader_extracts_ms1s(self):
        ms1s = read_mzml(self.args.mzml_files[0])
        with self.subTest('test that the right number of MS1s are found'):
            self.assertEqual(len(ms1s), 40)
        with self.subTest('test that all spectra are MS1s'):
            self.assertTrue(all([s[1].ms_level == 1 for s in ms1s]))

class processPSMtestSuite(base_test_classes.InitializedPSMsTestSuite):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.read_mzml = read_mzml
    
    def tearDown(self):
        super().tearDown()
        process_psm.__globals__['read_mzml'] = self.read_mzml
    
    def test_sinlge_psm_gets_needed_data(self):
        #initialize PSM list    
        PSM_list = [self.psms[0]]
        PSM_list[0].scan = 10
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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            results = process_psm('test1.mzML')
        with self.subTest('test that only one result is returned'):
            self.assertEqual(len(results), 1)
        with self.subTest('test that psm extracted enough peaks'):
            N_real_peaks = len([p for p in results[0].intensity if np.isfinite(p)])
            self.assertEqual(N_real_peaks, len(spectrum))
        
    def test_single_psm_ignores_noise(self):
        rng = np.random.default_rng(1)
        #initialize PSM list    
        PSM_list = [self.psms[0]]
        PSM_list[0].scan = 10
        #simulate unlabeled spectrum
        spectrum = isotopic_variants(PSM_list[0].formula, npeaks = 6, charge = PSM_list[0].charge)
        mz = [p.mz for p in spectrum]
        intensity = [p.intensity for p in spectrum]
        good_peaks = set(zip(mz, intensity, strict = True))
        N = 100
        mz += list(rng.uniform(0,100,N))
        intensity += list(rng.uniform(0,100,N))
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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            results = process_psm('test1.mzML')
        
        with self.subTest('test that the right number of peaks are extracted'):
            N_real_peaks = len([p for p in results[0].intensity if np.isfinite(p)])
            self.assertEqual(N_real_peaks, len(good_peaks))
        
        for peak in list(zip(results[0].mz, results[0].intensity, strict = True))[:len(good_peaks)]:
            closest = min(good_peaks, key = lambda x: abs(x[0] - peak[0]))
            with self.subTest('test mz close enough'):
                self.assertAlmostEqual(closest[0], peak[0])
            with self.subTest('test intensity close enough'):
                self.assertAlmostEqual(closest[1], peak[1])
    
    def test_multiple_psms_work(self):
        rng = np.random.default_rng(1)
        #initialize PSM list    
        PSM_list = [self.psms[0]]*5
        PSM_list[0].scan = 10
        #simulate unlabeled spectrum
        spectrum = isotopic_variants(PSM_list[0].formula, npeaks = 6, charge = PSM_list[0].charge)
        mz = [p.mz for p in spectrum]
        intensity = [p.intensity for p in spectrum]
        good_peaks = set(zip(mz, intensity, strict = True))
        N = 100
        mz += list(rng.uniform(0,100,N))
        intensity += list(rng.uniform(0,100,N))
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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            results = process_psm('test1.mzML')
        
        with self.subTest('test there are the right number of results'):
            self.assertEqual(len(results), 5)
        for result in results:
            with self.subTest('test each psm has the right number of peaks'):
                N_real_peaks = len([p for p in results[0].intensity if np.isfinite(p)])
                self.assertEqual(N_real_peaks, len(good_peaks))

class processSpectrumDataTestSuite(base_test_classes.InitializedPSMsTestSuite):
    def test_collection_of_psms(self):
        rng = np.random.default_rng(1)
        #initialize PSM list    
        PSM_list = [self.psms[0]]*5
        PSM_list[0].scan = 10
        #simulate unlabeled spectrum
        spectrum = isotopic_variants(PSM_list[0].formula, npeaks = 6, charge = PSM_list[0].charge)
        mz = [p.mz for p in spectrum]
        intensity = [p.intensity for p in spectrum]
        good_peaks = set(zip(mz, intensity, strict = True))
        N = 100
        mz += list(rng.uniform(0,100,N))
        intensity += list(rng.uniform(0,100,N))
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
            results = process_spectrum_data(self.args, PSM_list)
        
        with self.subTest('test there are the right number of results'):
            self.assertEqual(len(results), 5)
        for result in results:
            with self.subTest('test each psm has the right number of peaks'):
                N_real_peaks = len([p for p in results[0].intensity if np.isfinite(p)])
                self.assertEqual(N_real_peaks, len(good_peaks))

if __name__ == '__main__':
    unittest.main()
