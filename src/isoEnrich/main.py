#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 17:09:56 2024

@author: 4vt
"""

from isoEnrich.options import options
args = options()

from isoEnrich.parse_mzml import parse_PD, initialize_psms, process_spectrum_data
from isoEnrich.make_peptides import initialize_peptides

psm_data = parse_PD(args)
psms = initialize_psms(args, psm_data)
psms = process_spectrum_data(args, psms)
peptides = initialize_peptides(psms)
