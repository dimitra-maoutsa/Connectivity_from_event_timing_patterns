#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 02:25:09 2018

"otsu_method.py" : Calculates optimal mode seperating threshold for samples arising from a bimodal 
distribution that maximises the inter-mode variance according to Otsu method
[Nobuyuki Otsu (1979). "A threshold selection method from gray-level histograms".
 IEEE Trans. Sys., Man., Cyber. 9 (1): 62–66.]

Input
----------
>>> "X" : samples from bimodal distribution 

Output
----------
>>> "threshold": optimal threshold that seperates the 2 distributions according to Otsu method

Accompanying material to "Inferring network connectivity from event timing
patterns".


@author: dimitra
"""
import numpy as np

def otsu(X):
    N = X.size
    nbins = int(np.floor(0.5*N))
    counts,bis = np.histogram(X,nbins);
    p = counts / np.sum(counts);
    sigma_b = np.zeros(nbins)
    for t in range(nbins):
       q_L = np.sum(p[: t+1])
       q_H = np.sum(p[t+1: ])
       miu_L = np.sum(p[ : t] * np.arange(1,t+1).T) / q_L;
       miu_H = np.sum(p[t+1 : ] * np.arange(t + 1,nbins).T) / q_H;
       sigma_b[t] = q_L * q_H * (miu_L - miu_H)**2;
    
    
    threshold_otsu = np.argmax(sigma_b[:-1]);
    threshold = bis[threshold_otsu+1]
    return threshold