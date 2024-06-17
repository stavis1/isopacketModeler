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

from isopacketModeler.options import options


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
        
