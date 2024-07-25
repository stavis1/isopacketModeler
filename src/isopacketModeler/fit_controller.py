#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 11:29:53 2024

@author: 4vt
"""

from multiprocessing import Pool, Manager
import traceback

from isopacketModeler.data_generating_processes import BetabinomQuiescentMix, Betabinom, BinomQuiescentMix, Binom

data_generating_processes = {'BetabinomQuiescentMix':BetabinomQuiescentMix,
                             'Betabinom':Betabinom,
                             'BinomQuiescentMix':BinomQuiescentMix,
                             'Binom':Binom}

class peptide_fit_conroller():
    def __init__(self, args):
        self.DGPs = [data_generating_processes[p](args) for p in args.data_generating_processes]
        self.cores = args.cores

    def model_selection(self, peptide):
            peptide.canonical_fit = min(peptide.fit_results, key = lambda x: x.fit)

    def fit_all_DGPs(self, i):
        if event.is_set():
            return
        try:
            results = []
            for DGP in self.DGPs:
                results.append(DGP.fit(all_peptides[i]))
            return results
        except:
            traceback.print_exc()
            event.set()

    def fit_peptides(self, peptides):
        global all_peptides
        all_peptides = peptides
        
        with Manager() as manager:
            shared_event = manager.Event()
            def init_worker(shared_event):
                global event
                event = shared_event
            
            with Pool(processes = self.cores,
                      initializer=init_worker, 
                      initargs=(shared_event,)) as p:
                results = p.map(self.fit_all_DGPs, range(len(all_peptides)))
            
        for peptide, result in zip(peptides, results):
            peptide.fit_results = result
        
        for peptide in peptides:
            self.model_selection(peptide)

