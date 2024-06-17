#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 18:14:31 2024

@author: anon
"""

import unittest

import pandas as pd

from base_test_classes import ParsedOptionsTestSuite
from isopacketModeler.parse_mzml import parse_PSMs, initialize_psms, process_psm, process_spectrum_data

class parsePSMsTestSuite(ParsedOptionsTestSuite):
    def test_parses_PD_psms_file(self):
        psms = parse_PSMs(self.args)
        with self.subTest('Test that 6 psms were found'):
            self.assertEqual(psms.shape[0], 10)
        for psm in psms:
            with self.subTest('Test that PSM have the right number of columns'):
                self.assertEqual(psms.shape[1], 5)
    
    def test_parser_collects_arbitrary_metadata(self):
        self.args.PSM_headers += ["PSMs Workflow ID", "Tags"]
        psms = parse_PSMs(self.args)
        for psm in psms:
            with self.subTest('Test that PSM tuples are the right length'):
                self.assertEqual(psms.shape[1], 7)

class initializePSMsTestSuite(ParsedOptionsTestSuite):
    def test_label_file_PSMs_initialize_correctly(self):
        N = 3
        seq = ['TEST', 'A.TEST.A', 'TES[1.23]T']
        file = ['test1.mzML']*N
        scan = range(N)
        charge = [2]*N
        prots = ['TESTPROT']*N
        psm_data = pd.DataFrame({'seq':seq,
                                 'file':file,
                                 'scan':scan,
                                 'charge':charge,
                                 'proteins':prots})
        psms = initialize_psms(self.args, psm_data)
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
            self.assertListEqual([p.label for p in psms], ['C']*N)
    
    def test_control_PSMs_are_duplicated(self):
        N = 3
        seq = ['TEST', 'A.TEST.A', 'TES[1.23]T']
        file = ['test2.mzML']*N
        scan = range(N)
        charge = [2]*N
        prots = ['TESTPROT']*N
        psm_data = pd.DataFrame({'seq':seq,
                                 'file':file,
                                 'scan':scan,
                                 'charge':charge,
                                 'proteins':prots})
        psms = initialize_psms(self.args, psm_data)
        with self.subTest('Check the number of PSMS'):
            self.assertEqual(len(psms), N*2)
        with self.subTest('Check both labels are used'):
            self.assertListEqual([p.label for p in psms], ['N']*N + ['C']*N)
            
            
if __name__ == '__main__':
    unittest.main()




