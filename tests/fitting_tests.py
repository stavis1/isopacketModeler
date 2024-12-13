#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:33:50 2024

@author: 4vt
"""

import unittest

import numpy as np

from isopacketModeler.data_generating_processes import BetabinomQuiescentMix, Betabinom, BinomQuiescentMix, Binom
from isopacketModeler.fit_controller import peptide_fit_conroller, data_generating_processes
from isopacketModeler.data_objects import peptide
import base_test_classes

class BetabinomQuiescentMixTestSuite(base_test_classes.DataGeneratingProcessTestSuite):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.reasonable_params = np.array([0.5,4,3])
    
    def setUp(self):
        super().setUp()
        self.DGP = BetabinomQuiescentMix(self.args)

class BetabinomTestSuite(base_test_classes.DataGeneratingProcessTestSuite):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.reasonable_params = np.array([4,3])
    
    def setUp(self):
        super().setUp()
        self.DGP = Betabinom(self.args)

class BinomQuiescentMixTestSuite(base_test_classes.DataGeneratingProcessTestSuite):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.reasonable_params = np.array([0.5,0.5])
    
    def setUp(self):
        super().setUp()
        self.DGP = BinomQuiescentMix(self.args)

class BinomTestSuite(base_test_classes.DataGeneratingProcessTestSuite):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.reasonable_params = np.array([0.5])
    
    def setUp(self):
        super().setUp()
        self.DGP = Binom(self.args)

class FitControllerTestSuite(base_test_classes.InitializedPSMsTestSuite):
    def setUp(self):
        super().setUp()
        self.rng = np.random.default_rng(1)
        
    def make_peptide(self, DGP, reasonable_params):
        for psm in self.psms:
            psm.intensity = self.rng.uniform(1e5,1e7,len(psm.mz))
            ppm5 = (psm.mz[0]/1e6)*5
            psm.mz_err = self.rng.uniform(-ppm5,ppm5,len(psm.mz))
        pep = peptide(self.psms)
        expected_intensity = DGP.expected(pep, reasonable_params)
        expected_intensity = np.array([expected_intensity]*len(self.psms))
        pep.obs = expected_intensity
        return pep
    
    def test_fitting_multiple_peptides_works(self):
        reasonable_params = np.array([0.5,4,3])
        peptides = [self.make_peptide(BetabinomQuiescentMix(self.args), reasonable_params) for _ in range(10)]
        fit_controller = peptide_fit_conroller(self.args)
        peptides = fit_controller.fit_peptides(peptides)
        
        with self.subTest('are the correct number of peptides returned'):
            self.assertEqual(len(peptides), 10)
        
        for pep in peptides:
            with self.subTest('does each peptide have 4 results'):
                self.assertEqual(len(pep.fit_results), len(self.args.data_generating_processes))
            for dgp in pep.fit_results:
                with self.subTest('were all DGPs fit'):
                    self.assertIn(dgp.dgp_name, self.args.data_generating_processes)
            with self.subTest('was the correct DGP chosen'):
                self.assertEqual(pep.canonical_fit.dgp_name, 'BetabinomQuiescentMix')
            with self.subTest('were the fits good'):
                self.assertLess(pep.canonical_fit.fit, 0.01)
            for truth, fitted in zip(np.array([0.5,4,3]), pep.canonical_fit.params):
                with self.subTest('werer the correct parameters recovered'):
                    self.assertAlmostEqual(truth, fitted, delta = 0.001)
    
    def test_model_selection_recovers_model(self):
        DGPs = ['BetabinomQuiescentMix',
                'Betabinom',
                'BinomQuiescentMix',
                'Binom']
        params = [np.array([0.5,4,3]),
                  np.array([4,3]),
                  np.array([0.5,0.5]),
                  np.array([0.5])]
        
        class mock_event():
            def is_set(self):
                return False
        event = mock_event()
        fit_controller = peptide_fit_conroller(self.args)
        fit_controller.fit_all_DGPs.__globals__['event'] = event
        
        global all_peptides
        for DGP_name, param in zip(DGPs, params):
            DGP = data_generating_processes[DGP_name](self.args)
            pep = self.make_peptide(DGP, param)
            all_peptides = [pep]
            fit_controller.fit_all_DGPs.__globals__['all_peptides'] = all_peptides
            pep.fit_results = fit_controller.fit_all_DGPs(0)
            fit_controller.model_selection(pep)
            
            with self.subTest(f'was the correct DGP chosen for {DGP.name}'):
                self.assertEqual(pep.canonical_fit.dgp_name, DGP.name)
            with self.subTest(f'is the canonical fit good for {DGP.name}'):
                self.assertLess(pep.canonical_fit.fit, 0.01)
            with self.subTest(f'is the correct model fit good for {DGP.name}'):
                result = next(fit for fit in pep.fit_results if fit.dgp_name == DGP.name)
                self.assertLess(result.fit, 0.01)


if __name__ == '__main__':
    unittest.main()

