#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 17:09:56 2024

@author: 4vt
"""
from isopacketModeler.options import options
args = options()

from isopacketModeler.parse_mzml import parse_PSMs, initialize_psms, process_spectrum_data
from isopacketModeler.make_peptides import initialize_peptides
from isopacketModeler.classifier_tools import classifier
from isopacketModeler.fit_controller import peptide_fit_conroller
from isopacketModeler.report import make_report
from isopacketModeler.checkpoint import Checkpointer

checkpointer = Checkpointer(args)

#Collect data from mzML files.
if checkpointer.load_step < 1:
    psm_data = parse_PSMs(args)
    psms = initialize_psms(args, psm_data)
    psms = process_spectrum_data(args, psms)
    checkpointer.dump(psms, 1)
elif checkpointer.load_step == 1:
    psms = checkpointer.data

#Filter out unenriched PSMs.
if checkpointer.load_step < 2:
    bad_psms = []
    if args.do_psm_classification:
        psm_classifier = classifier(args)
        psm_data, psm_labels = psm_classifier.preprocess(psms)
        psm_classifier.fit(psm_data, psm_labels)
        psms, bad_psms = psm_classifier.winnow(psm_data, psms)
    
    #Collect PSMs into peptide objects.
    peptides = initialize_peptides(args, psms, bad_psms)
    checkpointer.dump(peptides, 2)
elif checkpointer.load_step == 2:
    peptides = checkpointer.data

#Fit models to peptide data.
fit_controller = peptide_fit_conroller(args)
peptides = fit_controller.fit_peptides(peptides)

#Export results.
make_report(args, peptides)

