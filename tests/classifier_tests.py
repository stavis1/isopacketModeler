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
import base_test_classes

class ClassifierTestSuite(base_test_classes.ParsedOptionsTestSuite):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.N = 10000
    
    def setUp(self):
        super().setUp()
        self.rng = np.random.default_rng(1)
        self.model = classifier_model(self.args)
    
    def make_false_data(self, rng):
        return rng.normal(0, 1, (256,2,self.N))

    def make_true_data(self, rng):
        return rng.normal(0.5, 1, (256,2,self.N))
        
    def test_fdr_control_works(self):
        def exp_fdr(x, m1, m2, s1, s2):
            return (1-norm.cdf(x, loc = m1, scale = s1))/norm.cdf(x, loc = m2, scale = s2)

        def err(x, m1, m2, s1, s2, target_FDR):
            return np.square(exp_fdr(x, m1, m2, s1, s2) - target_FDR)

        x_0 = 0
        m1, s1 = 0, 1
        m2, s2 = 2, 1
        target_FDR = self.model.fdr
        
        result = minimize(err, 
                          x0 = x_0,
                          args = (m1,m2,s1,s2,target_FDR))
        target_cutoff = result.x
        
        decoys = self.rng.normal(m1, s1, self.N)
        targets = np.concat((self.rng.normal(m1, s1, self.N), self.rng.normal(m2, s2, self.N)))
        self.model._set_cutoff(targets = targets, decoys = decoys)
        
        test_targets = np.concat((self.rng.normal(m1, s1, self.N), self.rng.normal(m2, s2, self.N)))
        ground_truth = np.array([0]*self.N + [1]*self.N)
        hits = ground_truth[test_targets > self.model.cutoff]
        test_FDR = np.sum(hits == 0)/len(hits)
        
        with self.subTest('test FDR is well controlled'):
            self.assertAlmostEqual(test_FDR, target_FDR, delta = 0.03)
        with self.subTest('test that the cutoff is close to where it should be'):
            self.assertAlmostEqual(self.model.cutoff, target_cutoff)

    def test_model_works_on_separable_data(self):
        ctrl = self.make_false_data(self.rng)
        obs = np.concatenate((self.make_false_data(self.rng),
                              self.make_true_data(self.rng)), axis = 2)
        data = np.concatenate((ctrl, obs), axis = 2)
        label = np.array([0]*self.N + [1]*(2*self.N))
        shuffle_idx = self.rng.choice(range(len(label)), len(label), replace = False)
        data = data[:,:,shuffle_idx]
        label = label[shuffle_idx]
        
        self.model.fit(data, label)
        
        new_obs = np.concatenate((self.make_false_data(self.rng),
                                  self.make_true_data(self.rng)), axis = 2)
        new_label = np.array([0]*self.N + [1]*self.N)
        calls = self.model._predict(new_obs)
        fdr = np.sum(np.logical_and(calls == 1, new_label == 0))/np.sum(calls)
        
        self.assertAlmostEqual(fdr, self.model.fdr, delta = 0.03)

    def test_model_excludes_inseparable_data(self):
        ctrl = self.make_false_data(self.rng)
        obs = np.concatenate((self.make_false_data(self.rng),
                              self.make_false_data(self.rng)), axis = 2)
        data = np.concatenate((ctrl, obs), axis = 2)
        label = np.array([0]*self.N + [1]*(2*self.N))
        shuffle_idx = self.rng.choice(range(len(label)), len(label), replace = False)
        data = data[:,:,shuffle_idx]
        label = label[shuffle_idx]
        
        self.model.fit(data, label)
        
        new_obs = np.concatenate((self.make_false_data(self.rng),
                                  self.make_false_data(self.rng)), axis = 2)
        calls = self.model._predict(new_obs)
        self.assertLess(np.sum(calls)/len(calls), 0.01)
    
    def test_preprocessor_works(self):
        class PSM:
            def __init__(self, rng, label):
                self.N = rng.integers(10,30)
                self.intensity = rng.uniform(0,1e6,self.N)
                self.mz_err = rng.uniform(-1,1,self.N)
                missing = rng.choice(range(self.N), 
                                     int(rng.uniform(0,0.9)*self.N), 
                                     replace = False)
                self.intensity[missing] = np.nan
                self.mz_err[missing] = np.nan
                self.mz = np.array(range(self.N))*rng.uniform(200,2000)
                self.label = label

        psms = [PSM(self.rng, '') for _ in range(500)] + [PSM(self.rng, 'C') for _ in range(500)]
        processed_data, processed_labels = self.model.preprocess(psms)
        
        with self.subTest('test that the right number of PSMs are preprocessed'):
            self.assertEqual(len(processed_labels), 1000)
        with self.subTest('test that the right number of PSMs are preprocessed'):
            self.assertEqual(len(processed_data), 1000)
        with self.subTest('ensure the data are shuffled'):
            self.assertAlmostEqual(np.sum(processed_labels[-500:])/1000, 0.5, delta = 0.2)

if __name__ == '__main__':
    unittest.main()





