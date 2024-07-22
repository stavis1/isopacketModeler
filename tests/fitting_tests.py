#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:33:50 2024

@author: 4vt
"""

import unittest

import numpy as np

from isopacketModeler.fitting_tools import BetabinomQuiescentMix
import base_test_classes

class BetabinomQuiescentMixTestSuite(base_test_classes.DataGeneratingProcessTestSuite):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.reasonable_params = np.array([0.5,4,3])
    
    def setUp(self):
        super().setUp()
        self.DGP = BetabinomQuiescentMix(self.args)

if __name__ == '__main__':
    unittest.main()

