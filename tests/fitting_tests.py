#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:33:50 2024

@author: 4vt
"""

import unittest

import numpy as np

from isopacketModeler.data_generating_processes import BetabinomQuiescentMix, Betabinom, BinomQuiescentMix, Binom
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



if __name__ == '__main__':
    unittest.main()

