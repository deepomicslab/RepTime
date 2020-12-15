import pandas as pd
import gzip
import sys
import time

def showTime(text):
    localtime = time.asctime( time.localtime(time.time()) )
    print("{0}: {1}\n".format(text, localtime), file=sys.stderr, flush=True)

def GC_bias(file_list):
    df = pd.DataFrame()
    for gc_file in file_list:
        df = df.append(pd.read_csv(gc_file,compression='gzip',sep="\t",
                                   usecols=[2,4],header=None, names=['RD','GC']))
    rd_avg = df.RD.mean()
    gc_avg = rd_avg / df.groupby('GC')['RD'].mean()
    gc_hash = gc_avg.to_dict()
    return(gc_hash)

def GC_norm(gc_file, gc_hash):
    #showTime('GC normlization start at')
    fh_GC = gzip.open(gc_file, 'rt')
    norm_file = gc_file.replace('ready','norm')
    fho = gzip.open(norm_file, 'wt')  
    for line in fh_GC:
        temp = line.strip().split("\t")
        gcnorm = int(temp[2]) * gc_hash[float(temp[4])]
        print("\t".join(temp+["%.2f"%gcnorm]), file = fho)
    #showTime('GC normlization start at')
        