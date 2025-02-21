#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:34:30 2024

@author: 4vt
"""

from scipy.optimize import basinhopping
from scipy.stats import betabinom, binom
import numpy as np

rng = np.random.default_rng(1)

class results():
    def __init__(self, dgp, peptide, minimizer_results):
        self.__dict__.update({k:v for k,v in minimizer_results.__dict__.items() if not k.startswith('__')})
        self.params = minimizer_results.x
        self.fit = minimizer_results.fun
        self.dgp_name = dgp.name
        self.fitted_dist = dgp.expected(peptide, minimizer_results.x)
        self.__dict__.update(dgp.extra_data(peptide, minimizer_results))

class DataGeneratingProcess():
    def __init__(self, args):
        self.name = NotImplemented
        self.bounds = NotImplemented
        #figure out how to pass initial guesses to get_x0
    
    def expected(self, peptide, params):
        raise NotImplementedError()
    
    def loss(self, params, peptide):
        exp = self.expected(peptide, params)
        #absolute error
        resids = np.abs(exp[np.newaxis, :] - peptide.obs)
        #winsorize at the 90th percentile
        q9 = np.nanquantile(resids, 0.9)
        resids[resids >= q9] = q9
        return np.nanmean(resids)
    
    def get_x0(self, peptide):
        raise NotImplementedError()
    
    def extra_data(self, peptide, result):
        return {}

    def fit(self, peptide):
        x_0 = self.get_x0(peptide)
        args = {'args':(peptide,), 
                # 'method':'SLSQP',#'Powell', 
                'bounds':self.bounds}
        result = basinhopping(self.loss, 
                              x_0,
                              T=2,
                              minimizer_kwargs = args,
                              take_step = self.take_step)
        return results(self, peptide, result)

class BBRandomDisplacementBounds(object):
    """random displacement with bounds for betabinomial models"""
    def __init__(self, bounds, stepsize=0.5):
        self.xmin = np.array([b[0] for b in bounds])
        self.xmax = np.array([b[1] for b in bounds])
        self.stepsize = stepsize

    def __call__(self, x):
        """take a random step but ensure the new position is within the bounds"""
        xnew = np.clip(x + rng.uniform(-self.stepsize, 
                                       self.stepsize, 
                                       np.shape(x)),
                       self.xmin,
                       self.xmax)
        #this boundry limits the mean enrichment to >= 0.05
        #the a parameter must be second to last and the b parameter must be last
        xnew[-2] = np.clip(xnew[-2], xnew[-1]/19, self.xmax[-2])
        return xnew

class BetabinomQuiescentMix(DataGeneratingProcess):
    def __init__(self, args):
        super().__init__(args)
        self.name = 'BetabinomQuiescentMix'
        self.bounds = [(0,1),
                       (1,100),
                       (1,100)]
        self.take_step = BBRandomDisplacementBounds(self.bounds)
    
    def extra_data(self, peptide, result):
        frac, a, b = result.x
        mean = a/(a+b)
        variance = (a*b)/(((a+b)**2)*(a+b+1))
        return {'mean':mean, 'variance':variance}
    
    def get_x0(self, peptide):
        return [0.5,4,3]
    
    def expected(self, peptide, params):
        label = betabinom.pmf(k = range(peptide.formula[peptide.label_elm]),
                              n = peptide.formula[peptide.label_elm],
                              a = params[1],
                              b = params[2])
        exp = np.convolve(label, peptide.background)
        exp = peptide.reshape(exp)
        exp = (peptide.unenriched*(1-params[0])) + (exp*params[0])
        exp = exp/np.nansum(exp)
        return exp

class Betabinom(DataGeneratingProcess):
    def __init__(self, args):
        super().__init__(args)
        self.name = 'Betabinom'
        self.bounds = [(1,100),
                       (1,100)]
        self.take_step = BBRandomDisplacementBounds(self.bounds)
    
    def extra_data(self, peptide, result):
        a, b = result.x
        mean = a/(a+b)
        variance = (a*b)/(((a+b)**2)*(a+b+1))
        return {'mean':mean, 'variance':variance}
    
    def get_x0(self, peptide):
        return [4,3]

    def expected(self, peptide, params):
        label = betabinom.pmf(k = range(peptide.formula[peptide.label_elm]),
                              n = peptide.formula[peptide.label_elm],
                              a = params[0],
                              b = params[1])
        exp = np.convolve(label, peptide.background)
        exp = peptide.reshape(exp)
        exp = exp/np.nansum(exp)
        return exp

class BinomRandomDisplacementBounds(object):
    """random displacement with bounds for betabinomial models"""
    def __init__(self, bounds, stepsize=0.5):
        self.xmin = np.array([b[0] for b in bounds])
        self.xmax = np.array([b[1] for b in bounds])
        self.stepsize = stepsize

    def __call__(self, x):
        """take a random step but ensure the new position is within the bounds"""
        xnew = np.clip(x + rng.uniform(-self.stepsize, 
                                       self.stepsize, 
                                       np.shape(x)),
                       self.xmin,
                       self.xmax)
        return xnew

class BinomQuiescentMix(DataGeneratingProcess):
    def __init__(self, args):
        super().__init__(args)
        self.name = 'BinomQuiescentMix'
        self.bounds = [(0,1),
                       (0.1,1)]
        self.take_step = BinomRandomDisplacementBounds(self.bounds)
    
    def get_x0(self, peptide):
        return [0.5,0.5]
    
    def expected(self, peptide, params):
        label = binom.pmf(k = range(peptide.formula[peptide.label_elm]),
                          n = peptide.formula[peptide.label_elm],
                          p = params[1])
        exp = np.convolve(label, peptide.background)
        exp = peptide.reshape(exp)
        exp = (peptide.unenriched*(1-params[0])) + (exp*params[0])
        exp = exp/np.nansum(exp)
        return exp

class Binom(DataGeneratingProcess):
    def __init__(self, args):
        super().__init__(args)
        self.name = 'Binom'
        self.bounds = [(0,1)]
        self.take_step = BinomRandomDisplacementBounds(self.bounds)
    
    def get_x0(self, peptide):
        return [0.5]
    
    def expected(self, peptide, params):
        label = binom.pmf(k = range(peptide.formula[peptide.label_elm]),
                          n = peptide.formula[peptide.label_elm],
                          p = params[0])
        exp = np.convolve(label, peptide.background)
        exp = peptide.reshape(exp)
        exp = exp/np.nansum(exp)
        return exp