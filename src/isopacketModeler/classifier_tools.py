#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 13:29:55 2024

@author: 4vt
"""
from copy import copy
from collections import defaultdict, Counter
import os
import logging

import numpy as np
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
logging.getLogger('tensorflow').setLevel(logging.ERROR)
os.environ["KMP_AFFINITY"] = "noverbose"
import tensorflow as tf
tf.autograph.set_verbosity(3)
from tensorflow.keras import layers, models
tf.config.set_visible_devices([], 'GPU')

class classifier():
    def __init__(self, args):
        self.args = args
        self.FDR = args.classifier_fdr
        self.cutoff = np.nan
        self.rng = np.random.default_rng(1)
        self.history = defaultdict(lambda : [])

    def _get_model(self):
        if 'model' in self.__dict__.keys():
            del self.model
        model = models.Sequential()
        model.add(layers.Conv2D(5, (3, 2), activation='relu', input_shape=(256, 2, 1), padding='same'))
        model.add(layers.Dropout(0.1))
        model.add(layers.MaxPooling2D((2, 1)))
        model.add(layers.Conv2D(5, (4, 5), activation='relu', padding='same'))
        model.add(layers.Dropout(0.1))
        model.add(layers.MaxPooling2D((2, 1)))
        model.add(layers.Conv2D(5, (4, 5), activation='relu', padding='same'))
        model.add(layers.Flatten())
        model.add(layers.Dropout(0.1))
        model.add(layers.Dense(50, activation='relu'))
        model.add(layers.Dropout(0.1))
        model.add(layers.Dense(50, activation='relu'))
        model.add(layers.Dense(1, activation = 'sigmoid'))
        
        model.compile(optimizer='adam',
                      loss=tf.keras.losses.MeanSquaredError(),
                      metrics=['accuracy', tf.keras.metrics.RecallAtPrecision(0.99)])
        return model

    def _fit_one_step(self, X, y, epochs = 7):
        self.model = self._get_model()
        history = self.model.fit(X, y.reshape((-1,1)), epochs=epochs)
        for key in history.history:
            self.history[key].extend(history.history[key])
        self.history['epochs'].append(epochs)
        return self

    def _set_cutoff(self, targets, decoys):
        ndecoy = len(decoys)
        ntarget = len(targets)
        decoy_scale = ntarget/ndecoy
        ŷ = list(decoys) + list(targets)
        y = [0]*ndecoy + [1]*ntarget
        observations = sorted(zip(ŷ,y), reverse = True)
        exp_fdr = (ndecoy*decoy_scale)/ntarget
        while exp_fdr > self.FDR:
            elm = observations.pop()
            if elm[1]:
                ntarget -= 1
            else:
                ndecoy -= 1
            exp_fdr = (ndecoy*decoy_scale)/ntarget if ntarget > 0 else 0
        self.cutoff = elm[0]
    
    def predict_proba(self, X):
        return self.model.predict(X)
    
    def predict(self, X):
        probs = self.predict_proba(X)[:,0]
        return probs > self.cutoff
    
    def winnow(self, X, psms):
        classes = self.predict(X)
        psms = [p for p,c in zip(psms, classes) if c == 1 and p.is_labeled]
        bad_psms = [p for p,c in zip(psms, classes) if c != 1 and p.is_labeled]
        self.args.logs.info(f'{len(psms)} PSMs have passed the classifier model filter.')
        return (psms, bad_psms)

    def preprocess(self, psms):
        data = []
        for psm in psms:
            intensity = copy(psm.intensity)
            intensity = intensity/np.nanmax(intensity)
            intensity[np.isnan(intensity)] = np.zeros(intensity.shape)[np.isnan(intensity)]
            mz_err = copy(psm.mz_err)
            mz_err = (mz_err/psm.mz)*1e5
            mz_err[np.isnan(mz_err)] = np.full(mz_err.shape,-1)[np.isnan(mz_err)]
            x = np.linspace(np.nanmin(psm.mz),np.nanmax(psm.mz),256)
            interp_i = np.interp(x, psm.mz, intensity)
            interp_i = interp_i/np.nansum(interp_i)
            interp_mz = np.interp(x, psm.mz, mz_err)
            interp_mz = interp_mz/np.nansum(interp_mz)
            data.append(np.concatenate((interp_i[np.newaxis,:,np.newaxis],interp_mz[np.newaxis,:,np.newaxis]), axis = 2))
        data = np.concatenate(data, axis = 0)
        labels = np.array([psm.is_labeled for psm in psms])
        return (data, labels)
    
    def _update_y(self, X, y_init, y_current):
        ŷ = self.predict_proba(X)[:,0]
        if Counter([round(ŷ_i, 4) for ŷ_i in ŷ]).most_common(1)[0][1] > len(ŷ)/10:
            return y_current
        ŷ =  np.array([ŷ_i if yinit_i else yinit_i for ŷ_i,yinit_i in zip(ŷ, y_init)])
        return ŷ
    
    def fit(self, X, y, niter = 7):
        y = copy(y)

        #split off sample of elements to estimate the FDR cutoff
        fit_idx = self.rng.choice(range(X.shape[0]), int(X.shape[0]*0.8), replace = False)
        fit_idx_set = set(fit_idx)
        cutoff_idx = np.array([i for i in range(X.shape[0]) if not i in fit_idx_set])
        X_cut = X[cutoff_idx, :, :]
        y_cut = y[cutoff_idx]
        X = X[fit_idx, :, :]
        y = y[fit_idx]
        y_init = copy(y)
        
        self.args.logs.debug(f'Fitting started. There are {X.shape[0]} elements in the training dataset and {X_cut.shape[0]} elements in the FDR control set.')
        #run the I-EM algorithm
        self._fit_one_step(X, y)
        for _ in range(niter):
            y = self._update_y(X, y_init, y)
            self._fit_one_step(X, y)
        
        #do FDR control
        y = self.predict_proba(X_cut)[:,0]
        self._set_cutoff(targets = y[y_cut == 1], decoys = y[y_cut == 0])
        self.args.logs.info('The classifier model has been fit.')
        return self


