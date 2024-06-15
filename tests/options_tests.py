#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 12:07:49 2024

@author: 4vt
"""
import os
import sys
import unittest
from shutil import rmtree
from isopacketModeler.options import options, InputError

class optionsTestSuite(unittest.TestCase):
    def setUp(self):
        self.argv_init = sys.argv
        self.wd_init = os.getcwd()
        sys.argv = [sys.argv[0], '--options', 'options.toml']
        self.args = options()
    
    def tearDown(self):
        sys.argv = self.argv_init
        rmtree(self.args.output_directory)
        os.chdir(self.wd_init)
        
    def test_overwrite_protection(self):
        self.args.overwrite = False
        with self.assertRaises(FileExistsError):
            self.args.handle_working_directory()
    
    def test_input_validation(self):
        required = ['working_directory',
                    'output_directory',
                    'overwrite',
                    'log_file',
                    'log_level',
                    'design_file',
                    'mzml_dir',
                    'psms',
                    'cores']
        for attr in required:
            tmp = self.args.__dict__[attr]
            del self.args.__dict__[attr]
            with self.subTest(msg = f'Testing validation of {attr}'):
                with self.assertRaises(InputError):
                    self.args.validate_inputs()
            setattr(self.args, attr, tmp)
    
    def test_mzml_identification(self):
        with self.subTest(msg = 'Testing mzml list length'):
            self.assertEqual(len(self.args.mzml_files), 2)
        with self.subTest(msg = 'Testing basename list length'):
            self.assertEqual(len(self.args.base_names), 2)
    
    def test_mzml_exist(self):
        for mzml in self.args.mzml_files:
            with self.subTest(msg = f'Testing existance of {mzml}'):
                self.assertTrue(os.path.exists(mzml))

if __name__ == '__main__':
    unittest.main()
