import numpy as np
import pandas as pd  
import gzip          

def process_window(win_data):
    #print win_data
    #sys.exit()
    weight = len([_ for _ in win_data if _[1] > 0])
    win = pd.DataFrame(win_data)
    return np.array(win[1]).mean(), 1*weight

           
def slidingWindow(dp_file, outfile, full_win, slide_win,cutoff = 0.9):
    win_start = 1
    win_next = slide_win
    mywin = "%d_%d" % (win_start, win_next)
    win_depth = {}
    win_len = {}
    win_depth[mywin] = 0#defaultdict(int)
    win_len[mywin] = 0#defaultdict(int)
    win_list = [mywin]
    
    fh_norm = gzip.open(dp_file,'rt')
    fho = open(outfile,'wt')
    for line in fh_norm:
        line = line.strip().split("\t")
        pos, dp = int(line[1]), float(line[5])
        if  win_start <= pos <= win_next:
            win_depth[mywin] += dp
            win_len[mywin] += 1
        else:
            while win_next < pos:
                win_start += slide_win
                win_next += slide_win
                mywin = "%d_%d" % (win_start, win_next)
                win_depth[mywin] = 0
                win_len[mywin] = 0
                win_list.append(mywin)
            win_depth[mywin] += dp
            win_len[mywin] += 1
    #print sliding window result        
    mynum = int(full_win / slide_win)
    i = 0
    tol_len = len(win_list) - mynum + 1
    while i < tol_len:
        mywin = win_list[i]
        mydepth = 0
        mylen = 0 
        for x in range(i, i+mynum):
            mydepth +=  win_depth[win_list[x]]
            mylen += win_len[win_list[x]]          
        weight = round(mylen / full_win, 4)
        if weight >= cutoff:
            mypos = mywin.split('_')[1]
            avgdepth = round(mydepth / mylen, 2)
            print ("%s\t%s\t%s" % (mypos, avgdepth, weight) ,file = fho,flush=True)
        i += 1
