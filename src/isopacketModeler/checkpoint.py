#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 13:25:01 2025

@author: 4vt
"""
import os
import sys
import dill
import glob

class Checkpointer:
    def __init__(self, options):
        self.logs = options.logs
        self.output_directory = options.output_directory
        self.stopping_point = options.stopping_point
        self.checkpoints = []
        for file in options.checkpoint_files:
            if '*' in file:
                self.checkpoints.extend(glob.glob(file))
            else:
                self.checkpoints.append(file)
        if self.checkpoints:
            self.load_step, self.data = self.load()
            self.logs.info(f'Loaded checkpoint for step {self.load_step}.')
        else:
            self.load_step = 0
            self.data = None
        
    def load(self):
        steps = set()
        data = []
        for checkpoint_file in self.checkpoints:
            self.logs.debug(f'Now loading checkpoint file {checkpoint_file}')
            with open(checkpoint_file, 'rb') as dillfile:
                step, data_tmp = dill.load(dillfile)
            data.extend(data_tmp)
            steps.add(step)
            if len(steps) > 1:
                self.logs.error('Checkpoints from multiple steps detected. Please only load checkpoints from one step.')
                sys.exit(1)
        return next(s for s in steps), data
    
    def dump(self, data, step):
        checkpoint_file = f'{self.output_directory}checkpoint_step{step}_{str(hash(str(data))%9999) + str(os.getpid())}.dill'
        with open(checkpoint_file, 'wb') as dillfile:
                dill.dump((step, data), dillfile)
        self.logs.info(f'Saved checkpoint for step {step} in file {checkpoint_file}.')
        
        if self.stopping_point == step:
            self.logs.info(f'Now stopping at checkpoint step {step}.')
            sys.exit(0)
    
    