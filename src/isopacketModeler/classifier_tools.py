#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 13:29:55 2024

@author: 4vt
"""
from copy import copy
from collections import defaultdict
import os
import re

import dill
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow.keras import layers, models

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

def preprocess(intensity, mz_err, mz):
    intensity = copy(intensity)
    intensity = intensity/np.nanmax(intensity)
    intensity[np.isnan(intensity)] = np.zeros(intensity.shape)[np.isnan(intensity)]
    mz_err = copy(mz_err)
    mz_err = (mz_err/mz)*1e5
    mz_err[np.isnan(mz_err)] = np.full(mz_err.shape,-1)[np.isnan(mz_err)]
    x = np.linspace(min(mz),max(mz),256)
    interp_i = np.interp(x, mz, intensity)
    interp_mz = np.interp(x, mz, mz_err)
    return np.concatenate((interp_i[:,np.newaxis],interp_mz[:,np.newaxis]), axis = 1)

class classifier_model:
    def __init__(self):
        model = models.Sequential()
        model.add(layers.Conv2D(32, (3, 2), activation='relu', input_shape=(256, 2, 1), padding='same'))
        model.add(layers.Dropout(0.1))
        model.add(layers.MaxPooling2D((2, 1)))
        model.add(layers.Conv2D(64, (3, 2), activation='relu', padding='same'))
        model.add(layers.Dropout(0.1))
        model.add(layers.MaxPooling2D((2, 1)))
        model.add(layers.Conv2D(64, (3, 2), activation='relu', padding='same'))
        model.add(layers.Flatten())
        model.add(layers.Dropout(0.1))
        model.add(layers.Dense(64, activation='relu'))
        model.add(layers.Dense(1, activation = 'sigmoid'))

        model.compile(optimizer='adam',
                      loss=tf.keras.losses.BinaryCrossentropy(),
                      metrics=['accuracy', tf.keras.metrics.RecallAtPrecision(0.99)])
        self.model = model
        self.history = defaultdict(lambda : [])
        
    def fit(self, features, labels, path = '', split = 1, epochs = 5):
        if split < 1:
            split_idx = int(len(labels)*split)
            train_features = features[:split_idx]
            train_labels = labels[:split_idx]
            test_features = features[split_idx:]
            test_labels = labels[split_idx:]
            history = self.model.fit(train_features, train_labels, epochs=epochs, 
                                validation_data=(test_features, test_labels))
        else:
            history = self.model.fit(features, labels, epochs=epochs)
        
        for key in history.history:
            self.history[key].extend(history.history[key])
        self.history['epochs'].append(epochs)
        if path:
            if not os.path.exists(path):
                os.mkdir(path)
            self.model.save_weights(os.path.abspath(path)+'/model.ckpt')
            with open(os.path.abspath(path)+'/history.dill','wb') as dillfile:
                dill.dump(self.history, dillfile)
    
    def load(self, path):
        self.model.load_weights(os.path.abspath(path)+'/model.ckpt').expect_partial()
        with open(os.path.abspath(path)+'/history.dill','rb') as dillfile:
            self.history = dill.load(dillfile)
    
    def evaluate(self, features, labels):
        eval_metrics = self.model.evaluate(features, labels, return_dict = True)
        print('\n'.join([f'{k}: {v}' for k,v in eval_metrics.items()]))
        
        preds, calls = self.predict(features)
        t = [p for p,l in zip(preds, labels) if l]
        f = [p for p,l in zip(preds, labels) if not l]
        fig, ax = plt.subplots()
        bins = np.linspace(0,1,100)
        ax.hist(t, bins = bins, color = 'k', alpha = 0.5, label = 'Labeled')
        ax.hist(f, bins = bins, color = 'r', alpha = 0.5, label = 'Control')
        ax.legend()
        plt.show()
    
    def show_history(self):
        n_epochs = np.sum(self.history['epochs'])
        keys = list(self.history.keys())
        for key in keys:
            match = re.search(r'(.+)_\d+\Z', key)
            if match is not None and match.group(1) in keys:
                self.history[match.group(1)] = self.history[match.group(1)] + self.history[key] 
                del self.history[key]
        
        for key in [k for k in self.history.keys() if k != 'epochs']:
            fig, ax = plt.subplots()
            ax.plot(range(n_epochs), self.history[key], '.-k')
            if 'loss' in key:
                ax.set_yscale('log')
            ax.set_ylabel(key)
            ax.set_xlabel('epochs')
            ax.set_xlim(0,n_epochs-1)
            if len(self.history['epochs']) > 1:
                ylim = ax.get_ylim()
                for i in range(1,len(self.history['epochs'])):
                    boundry = np.sum(self.history['epochs'][:i])
                    lines = ax.plot([boundry-0.5]*2,ylim, '--r', linewidth = 0.5)
                lines[0].set_label('Training Rounds')
                ax.set_ylim(ylim)
                ax.legend()
            plt.show()
    
    def _calc_cutoff(self, targs, decoys):
        n_decoy = len(decoys)
        n_preds = len(targs)
        ratio = n_decoy/n_preds
        targ = 0.01
        all_preds = sorted([(p,True) for p in targs] + [(p,False) for p in decoys])
        seen = 0
        seen_decoy = 0
        for pred in all_preds:
            if pred[1]:
                seen += 1
            else:
                seen_decoy +=1
            fdr = ((n_decoy - seen_decoy)/ratio)/(n_preds - seen)
            if fdr < targ:
                cutoff = pred[0]
                break
        return cutoff

    def predict(self, features):
        preds = [p[0] for p in self.model.predict(features)] 
        cutoff = 0.8
        calls = [p > cutoff for p in preds]
        return (preds, calls)

