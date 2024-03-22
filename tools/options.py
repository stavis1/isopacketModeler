#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 15:03:53 2024

@author: 4vt
"""

from argparse import ArgumentParser
import tomllib
import os

class options:
    def __init__(self):
        parser = ArgumentParser()
        parser.add_argument('-o', '--options', action = 'store', required = True,
                            help = 'Path to options file.')
        args = parser.parse_args()
        
        with open(args.options,'rb') as toml:
            options = tomllib.load(toml)
        self.__dict__.update(options)
        
        required = ['working_directory',
                   'design_file',
                   'mzml_dir',
                   'psms']
        problems = [r for r in required if not r in self.__dict__.keys()]
        if problems:
            raise Exception('Required settings not found in options file:\n' + \
                            '\n'.join(problems))

        os.chdir(self.working_directory)
        self.find_mzml()
        self.parse_design()
        if not 'cores' in self.__dict__.keys():
            self.cores = 1
        
        self.dump()
    
    @classmethod
    def alt_init(cls):
        parser = ArgumentParser()
        parser.add_argument('-o', '--options', action = 'store', required = True,
                            help = 'Path to options file.')
        args = parser.parse_args()
        
        with open(args.options,'rb') as toml:
            options = tomllib.load(toml)
        
        os.chdir(options["working_directory"])
        if os.path.exists('options.dill'):
            import dill
            with open('options.dill','rb') as dillfile:
                return dill.load(dillfile)
        else:
            return cls()

    def find_mzml(self):
        mzml_files = [f for f in os.listdir(self.mzml_dir) if f.lower().endswith('.mzml')]
        self.mzml_files = [os.path.join(self.mzml_dir, f) for f in mzml_files]
        self.base_names = set(f[:-5] for f in mzml_files)
    
    def parse_design(self):
        import pandas as pd

        design = pd.read_csv(self.design_file, sep = '\t')
        design = design.fillna({'label':''})
        design['control'] = [not l for l in design['label']]
        self.design = design
        
        extra_design = [f for f in design['file'] if f not in self.base_names]
        if extra_design:
            print('Files specified in the design file without corresponding .mzML files:\n')
            print('\n'.join(extra_design))
        design_files = set(design['file'])
        extra_mzml = [f for f in self.base_names if f not in design_files]
        if extra_mzml:
            print('Files in the mzML directory without corresponding design entry files:\n')
            print('\n'.join(extra_mzml))
        
    def dump(self):
        import dill
        with open('options.dill', 'wb') as dillfile:
            dill.dump(self, dillfile)

