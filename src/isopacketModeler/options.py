#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 15:03:53 2024

@author: 4vt
"""

from argparse import ArgumentParser
import tomllib
import os
import logging

class InputError(Exception):
    pass

class options:
    def __init__(self):
        self.parse_args()
        self.handle_working_directory()
        self.logger_init()
        self.validate_inputs()        
        self.find_mzml()
        self.parse_design()
        if self.cores < 1:
            if 'SLURM_CPUS_PER_TASK' in os.environ.keys():
                self.cores = int(os.environ['SLURM_CPUS_PER_TASK'])
            else:
                self.cores = os.cpu_count()
        self.logs.debug(f'Now using {self.cores} cores.')
        self.parse_AA_formulae()
    
    def parse_args(self):
        parser = ArgumentParser()
        subparsers = parser.add_subparsers(help='Either specify options as a text file or as command line arguments.')
        file_parser = subparsers.add_parser('file')
        file_parser.add_argument('-o', '--options', action = 'store', required = True,
                                 help = 'Path to options file.')

        cmd_parser = subparsers.add_parser('cmd')
        cmd_parser.add_argument('--working_directory', action = 'store', required = True,
                                help = 'All paths should be relative to this directory')
        cmd_parser.add_argument('--output_directory', action = 'store', required = True,
                                help = 'Directory to store results')
        cmd_parser.add_argument('--design_file', action = 'store', required = True,
                                help = 'A TSV file defining experimental design; see README for details')
        cmd_parser.add_argument('--mzml_dir', action = 'store', required = True,
                                help = 'The top level directory under which all mzML files can be found')
        cmd_parser.add_argument('--psms', action = 'append', required = True,
                                help = 'Use once per file for all PSM files')
        cmd_parser.add_argument('--psm_headers', action = 'store', required = True,
                                help = 'A comma separated list of column names in the PSM files: e.g. sequence,mzml,scan#,charge,proteins')
        cmd_parser.add_argument('--AA_formulae', action = 'store', required = True,
                                help = 'The amino acid chemical formula file; see README for details')
        cmd_parser.add_argument('--cores', action = 'store', required = False, default = 0, type = int,
                                help = 'The maximum number of cores to use')
        cmd_parser.add_argument('--classifier_fdr', action = 'store', required = False, default = 0.05, type = float,
                                help = 'The false discovery rate target for the PSM classifier')
        cmd_parser.add_argument('--data_generating_processes', action = 'append', required = True, choices = ['BetabinomQuiescentMix',
			                                                                                                  'Betabinom',
   		                                                                                                      'BinomQuiescentMix',
                                                                                                              'Binom'],
                                help = 'Use once for each model to fit to peptide data')
        cmd_parser.add_argument('--max_peptide_err', action = 'store', required = False, type = float, default = 0.015,
                                help = 'Peptides with model fit error above this value will not be reported')
        cmd_parser.add_argument('--do_PSM_classification', action = 'store_true', required=False, default=False,
                                help = 'Whether to do a preliminary classification of isotope enrichment')
        cmd_parser.add_argument('--checkpoint_files', action = 'append', required = False,
                                help = 'Use once per checkpoint file, all checkpoints must be at the same step')
        cmd_parser.add_argument('--stopping_point', action = 'store', required = False, default = False, type = int, choices = [1,2],
                                help = 'What step to stop at if you wish to stop early')
        args = parser.parse_args()
        
        if hasattr(args, 'options'):
            with open(args.options,'rb') as toml:
                options = tomllib.load(toml)
            self.__dict__.update(options)
        else:
            self.__dict__.update({k:v for k,v in args.__dict__.items() if v is not None and not k.startswith('_')})
    
    def handle_working_directory(self):
        os.chdir(self.working_directory)
        if not os.path.exists(self.output_directory):
            os.mkdir(self.output_directory)      
        else:
            if not self.overwrite:
                raise FileExistsError('An output directory with this name already exists and overwrite is false.')
    
    def logger_init(self):
        self.logs = logging.getLogger('isopacketModeler')
        self.logs.setLevel(10)
        formatter = formatter = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s')

        logfile = logging.FileHandler(self.log_file)
        logfile.setLevel(10)
        logfile.setFormatter(formatter)
        self.logs.addHandler(logfile)
        
        logstream = logging.StreamHandler()
        logstream.setLevel(self.log_level)
        logstream.setFormatter(formatter)
        self.logs.addHandler(logstream)
    
    def validate_inputs(self):
        #check if toml has all necessary information
        required = ['working_directory',
                    'output_directory',
                    'overwrite',
                    'log_file',
                    'log_level',
                    'design_file',
                    'mzml_dir',
                    'psms',
                    'cores',
                    'AA_formulae',
                    'PSM_headers',
                    'data_generating_processes',
                    'max_peptide_err',
                    'do_PSM_classification']
        problems = [r for r in required if not r in self.__dict__.keys()]
        if problems:
            msg = 'Required settings not found in options file:\n' + '\n'.join(problems)
            self.logs.error(msg)
            raise InputError()
        if type(self.PSM_headers) == str:
            self.PSM_headers = self.PSM_headers.split(',')
            
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
            message = 'Files specified in the design file without corresponding .mzML files:\n'
            message += '\n'.join(extra_design)
            self.logs.warning(message)

        design_files = set(design['file'])
        extra_mzml = [f for f in self.base_names if f not in design_files]
        if extra_mzml:
            message = 'Files in the mzML directory without corresponding design entry files. These will be ignored:\n'
            message += '\n'.join(extra_mzml)
            self.logs.warning(message)
        self.mzml_files = [f for f in self.mzml_files if os.path.basename(f)[:-5] in design_files]
        self.base_names = [f for f in self.base_names if f in design_files]
        
    def parse_AA_formulae(self):
        with open(self.AA_formulae, 'r') as tsv:
            cols = tsv.readline().strip().split('\t')
            self.AA_formulae = {l.split('\t')[0]:{e:int(c) for e,c in zip(cols[1:],l.strip().split('\t')[1:], strict = True)} for l in tsv}
