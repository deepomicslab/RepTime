# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 13:52:36 2020

@author: WANG
"""

import sys
import numpy as np
import pandas as pd

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return average, np.sqrt(variance)

def getMeanSD(file_list, cent=0.8, wighted=False):
    df = pd.DataFrame()
    for win_file in file_list:
        df = df.append(pd.read_csv(win_file,sep="\t",header=None,usecols=[1,2]))
    wins = np.array(df)
    wins = wins[wins[:,0].argsort()]
    center_wins = wins[int(((1-cent)/2)*len(wins)): int((cent + (1-cent)/2)*len(wins))]  # not exact, since we are not considering weight
    if wighted:
        avg_trimmed, sd_trimmed = weighted_avg_and_std(center_wins[:,0], center_wins[:,1])
    else:
        avg_trimmed = np.mean(center_wins[:,0])
        sd_trimmed = np.std(center_wins[:,0])
    return(avg_trimmed, sd_trimmed)

def filterwindow(win_file, mean, sd, weight_cutoff):
    upper = mean + 3 * sd
    lower = mean - 3 * sd
    fh = open(win_file)
    outfile = win_file.replace('.raw', '')
    fho = open(outfile, 'wt')
    for l in fh:
        sp = l.strip().split()
        dp = float(sp[1])
        weight = float(sp[2])
        if lower <= dp <= upper and weight > weight_cutoff:
            #sys.stdout.write('{0}\t{1}\t{2}\n'.format(sp[0], (dp - mean)/sd, sp[2])) 
            #rint('{0}\t{1}\t{2}\t{3}'.format(sp[0], (dp - mean)/sd, sp[2], dp),file = fho)
            print('{0}\t{1}\t{2}'.format(sp[0], dp, sp[2]),file = fho)