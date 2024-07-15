#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 15:39:24 2024

@author: anon
"""

import unittest

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

from isopacketModeler.classifier_tools import classifier_model
from base_test_classes import ParsedOptionsTestSuite

class FDRcontrolTestSuite(ParsedOptionsTestSuite):
    def setUp(self):
        super().setUp()
        self.rng = np.random.default_rng(1)
    
    def test_fdr_cutoff_is_close(self):
        def exp_fdr(x, m1, m2, s1, s2):
            return (1-norm.cdf(x, loc = m1, scale = s1))/norm.cdf(x, loc = m2, scale = s2)

        def err(x, m1, m2, s1, s2, target_FDR):
            return np.square(exp_fdr(x, m1, m2, s1, s2) - target_FDR)
        
        model = classifier_model(self.args)

        x_0 = 0
        m1, s1 = 0, 1
        m2, s2 = 2, 1
        target_FDR = model.fdr
        
        result = minimize(err, 
                          x0 = x_0,
                          args = (m1,m2,s1,s2,target_FDR))
        target_cutoff = result.x
        
        N = 10000
        decoys = self.rng.normal(m1, s1, N)
        targets = np.concat((self.rng.normal(m1, s1, N), self.rng.normal(m2, s2, N)))
        model._set_cutoff(targets = targets, decoys = decoys)
        
        test_targets = np.concat((self.rng.normal(m1, s1, N), self.rng.normal(m2, s2, N)))
        ground_truth = np.array([0]*N + [1]*N)
        hits = ground_truth[test_targets > model.cutoff]
        test_FDR = np.sum(hits == 0)/len(hits)
        
        self.assertAlmostEqual(test_FDR, target_FDR, delta = 0.01)






