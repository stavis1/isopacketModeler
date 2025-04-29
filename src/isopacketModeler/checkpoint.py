#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 13:25:01 2025

@author: 4vt
"""
import os
import sys
import dill

class Checkpointer:
    def __init__(self, options):
        self.logs = options.logs
        self.output_directory = options.output_directory
        self.stopping_point = options.stopping_point
        self.checkpoints = options.checkpoint_files
        if self.checkpoints:
            self.load_step, self.data = self.load()
        else:
            self.load_step = 0
            self.data = None
        
    def load(self):
        steps = set()
        data = []
        for checkpoint_file in self.checkpoints:
            with open(checkpoint_file, 'rb') as dillfile:
                step, data_tmp = dill.load(dillfile)
            data.extend(data_tmp)
            steps.add(step)
            if len(steps) > 1:
                self.logs.error('Checkpoints from multiple steps detected. Please only load checkpoints from one step.')
                sys.exit(1)
        return next(s for s in steps), data
    
    def dump(self, data, step):
        with open(f'{self.output_directory}checkpoint_step{step}_{os.getpid()}.dill', 'wb') as dillfile:
            dill.dump((step, data), dillfile)
        
        if self.stopping_point == step:
            self.logs.error(f'Now stopping at checkpoint step {step}.')
            sys.exit(0)
    
    