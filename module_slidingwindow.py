# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 12:57:06 2020

@author: WANG
"""

import sys
import numpy as np
import pandas as pd            

def process_window(win_data):
    #print win_data
    #sys.exit()
    weight = len([_ for _ in win_data if _[1] > 0])
    win = pd.DataFrame(win_data)
    return np.array(win[1]).mean(), 1*weight

def slidingWindow(dp_file, full_win, slide_win,cutoff = 0.9):
    win_start = 0
    win_end = full_win
    allwin = []
    print_win = []
    outfile = dp_file.replace('gc.norm.gz', '10k.win.raw')
    fho = open(outfile,'wt')
    data = pd.read_csv(dp_file, sep='\t', header=None, compression='gzip',
                       usecols=[1, 5], names=['Pos', 'Rd'])
    for pos, dp in zip(data['Pos'], data['Rd']):
        allwin.append([pos,dp])
        if pos - win_end > slide_win:
            win_start = pos - 1
            win_end = pos + full_win - 1
        if pos >= win_end:
            for i in range(len(allwin)):
                if allwin[i][0] > win_start:
                    del allwin[:i]
                    break
            m, w = process_window(allwin)
            win_weight = w/full_win
            if(win_weight >= cutoff):
                base_pos = (allwin[0][0] + allwin[-1][0]) / 2
                print_win = [ str(x) for x in [base_pos, m, win_weight] ]
                print("\t".join(print_win), file=fho)
            win_start += slide_win
            win_end += slide_win