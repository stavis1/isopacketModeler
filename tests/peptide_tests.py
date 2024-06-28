#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 10:55:33 2024

@author: 4vt
"""

import base_test_classes
from isopacketModeler.make_peptides import fingerprint, initialize_peptides

class MakePepetidesTestSuite(base_test_classes.ProcessedPSMsTestSuite):
    def test_nonsimilar_PSMs_arent_merged(self):
        peptides = initialize_peptides(self.psms)
        with self.subTest('check the right number of peptides is generated'):
            self.assertEqual(len(peptides), len(self.psms))
        with self.subTest('check that all psms were included as peptides'):
            self.assertCountEqual([p.raw_sequence for p in self.psms], [p.raw_sequence for p in peptides])
    
    def test_fingerprint_identifies_unique_psms(self):
        psms = list(self.psms)*5
        fingerprints = [fingerprint(psm) for psm in psms]
        self.assertEqual(len(set(fingerprints)), len(self.psms))

    def test_similar_PSMs_are_merged(self):
        psms = list(self.psms)*5
        peptides = initialize_peptides(psms)
        with self.subTest('check the right number of peptides is generated'):
            self.assertEqual(len(peptides), len(self.psms))
        with self.subTest('check that all psms were included as peptides'):
            self.assertCountEqual([p.raw_sequence for p in self.psms], [p.raw_sequence for p in peptides])
    
