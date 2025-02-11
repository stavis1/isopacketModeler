#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 14:08:53 2024

@author: 4vt
"""

import dill
import pandas as pd

def make_report(args, peptides):
    #serialize peptide objects
    with open(f'{args.output_directory}peptides.dill', 'wb') as dillfile:
        dill.dump(peptides, dillfile)
    args.logs.debug(f'{len(peptides)} peptides have been saved.')
    
    #make fitted peptide table
    report = pd.DataFrame([p.report() for p in peptides])
    report.to_csv(f'{args.output_directory}peptides.tsv', sep = '\t', index = False)
    
