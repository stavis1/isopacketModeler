#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 18:14:31 2024

@author: anon
"""

import unittest

from base_test_classes import ParsedOptionsTestSuite
from isopacketModeler.parse_mzml import parse_PSMs, initialize_psms, process_psm, process_spectrum_data

class parsePSMsTestSuite(ParsedOptionsTestSuite):
    def test_parses_PD_psms_file(self):
        psms = list(parse_PSMs(self.args))
        with self.subTest('Test that 6 psms were found'):
            self.assertEqual(len(psms), 10)
        for psm in psms:
            with self.subTest('Test that PSM tuples are the right length'):
                self.assertEqual(len(psm), 5)
    
    def test_parser_collects_arbitrary_metadata(self):
        self.args.PSM_headers += ["PSMs Workflow ID", "Tags"]
        psms = list(parse_PSMs(self.args))
        for psm in psms:
            with self.subTest('Test that PSM tuples are the right length'):
                self.assertEqual(len(psm), 7)

if __name__ == '__main__':
    unittest.main()




