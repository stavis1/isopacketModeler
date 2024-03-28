#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 17:52:16 2024

@author: 4vt
"""
import sys
import os
import unittest
from options import options


class TestSuite(unittest.TestCase):
    def setUp(self):
        if '__file__' in globals():
            self.test_dir = os.path.dirname(os.path.abspath(__file__))
        else: 
            self.test_dir = '/home/4vt/Documents/data/SLT02_MagdiProject/software/enrichment_pipeline/tests/'
        sys.argv = [sys.argv[0],
                    '--options',
                    self.test_dir + 'options.toml']
        self.args = options()
    
    def tearDown(self):
        dillfiles = os.listdir(self.test_dir)
        for file in dillfiles:
            os.remove(self.test_dir + file)
        return
    
        
