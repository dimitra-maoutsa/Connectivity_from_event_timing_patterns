#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 18:54:29 2015

"calculate_AUCs.py" calculates the AUC score for identifying excitatory connections, 
inhibitory connections and weighted AUC score for identifying both E/I synapses

Input
--------
>>>> "estimates" : array with gradient estimates resulting from the regression problem  Eq. (9) in main article
>>>> "Jin"       : array with 1 designating existence of inhibitory synapses, and -1 absence of those
>>>> "Jex"       : array with 1 designating existence of excitatory synapses, and -1 absence of those 

Output
--------
>>>> "roc_aucin" : AUC score for inhibitory connections
>>>> "roc_aucex" : AUC score for excitatory connections
>>>> "wauc"      : weighted AUC score 
@author: Dimi
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc



def calculate_AUCs(estimates, Jin, Jex):    
    
    fpr, tpr, threx = roc_curve(Jex, estimates, pos_label=1)
    where_are_NaNs = np.isnan(fpr)
    fpr[where_are_NaNs] = 0
    where_are_NaNs = np.isnan(tpr)
    tpr[where_are_NaNs] = 0
    roc_aucex = auc(fpr, tpr)  
    fpr, tpr, thrin = roc_curve(Jin, estimates , pos_label=1)   
    
    where_are_NaNs = np.isnan(fpr)
    fpr[where_are_NaNs] = 0
    where_are_NaNs = np.isnan(tpr)
    tpr[where_are_NaNs] = 0
    roc_aucin = auc(fpr, tpr)
    
    
    countinh = len(list(filter(lambda a: a<0 , Jin))) #numb. of inhibitory synapses in ground truth
    countex = len(list(filter(lambda a: a>0 , Jex)))  #numb. of inhibitory synapses in ground truth
    
    P = countinh + countex
    #weighted AUC score according to the prevalence of the two types of interactions in ground truth (Jin,Jex)
    wauc = roc_aucin * (countinh/P) + roc_aucex * (countex/P)  
    return (roc_aucin, roc_aucex, wauc)
