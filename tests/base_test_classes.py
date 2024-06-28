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

import pandas as pd

from isopacketModeler.options import options
from isopacketModeler.parse_mzml import initialize_psms


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
        
