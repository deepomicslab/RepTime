# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 16:59:16 2020

@author: WANG
"""

import argparse
import os
import sys
import re
from multiprocessing import Process
from multiprocessing import Pool
import gzip
import time

parser = argparse.ArgumentParser(description="Process args")
parser.add_argument('-i', '--input', type = str, action='store', 
                    help='input bam file ')
parser.add_argument("-c", '--chr', type = str, action='store', default = 'chr1-22', 
                    help="chr for analysis.For example, chr1 or chr1-3,chr5,chr14 or 1-4,7,8-10,22")
parser.add_argument("-r", '--ref', type = str, action='store',
                    help="faidx indexed reference sequence file for input bam")
parser.add_argument("-o", '--outdir', type = str, action='store', default = './result', 
                    help="the directory of output result")
parser.add_argument("-p", '--process', type = int, action='store', default = 4, 
                    help="the number of process")
parser.add_argument("-s", '--step', type = str, action='store', default = '1,2,3,4', 
                    help="the step of repTime.1:mpileup; 2:GC; 3:sliding window; 4:filter; 5:smooth")
args = parser.parse_args()

binpath = sys.path[0]
cwd = os.getcwd()
if not args.outdir.startswith('/'):
    args.outdir = cwd + '/' + args.outdir
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

def showTime(text):
    localtime = time.asctime( time.localtime(time.time()) )
    print("{0}: {1}\n".format(text, localtime), file=sys.stderr, flush=True)
    
def chrList(chr):
    chr_list = []
    chr_type = re.match(r'(.*?)\d',chr).group(1)
    chr_string = chr.strip().replace(chr_type,'').split(',')
    for item in chr_string:
        if '-' in item:
            a, b = [ int(x) for x in item.split('-') ]
            my_chr = list(range(a, b+1))
        else:
            my_chr = [int(item)]
        chr_list += my_chr
    chrid_list = [ chr_type+str(x) for x in chr_list ]
    chrout_list = [ 'chr'+str(x) for x in chr_list ]
    return(chrid_list, chrout_list)

def mpileup(bam, ref, chr, outdir, samtools = binpath+'/samtools'):
    mpileup_file = "{0}/{1}.mpileup.txt.gz".format(outdir, chr)
    cmd = "{0} mpileup -f {1} -r {2} {3} |gzip > {4}".format(samtools, ref, chr, bam, mpileup_file)
    #print(cmd)
    #os.system(cmd)
    return(mpileup_file)
    

#def runRepTime(bam, chr, outdir, ref):
if __name__  == '__main__':
    bam = args.input
    chr = args.chr
    outdir = args.outdir
    ref = args.ref
    chr_list = []
    chr_list, out_list = chrList(chr)
    wing = 200
    unique = False
    #1# samtools mpileup
    if('1' in args.step):
        import module_mpileupGC 
        mpileup_pool = Pool(args.process)
        for chr_i in chr_list:
            ###multi process
            #mpileup(bam, ref, chr_i, outdir)
            mpileup_pool.apply_async(module_mpileupGC.runmpileupGC, args=(bam, ref, chr_i, outdir,wing,unique,))
            #runmpileupGC(bam, ref, chr, outdir, wing)
        mpileup_pool.close()
        mpileup_pool.join()
   
    #2# GC content normlization
    if('2' in args.step):
        import module_GCnorm
        ##2.1 hash table of GC content##
        showTime('GC correction start at')
        gc_file_list = [ os.path.join(outdir, x+'.gc.ready.gz') for x in out_list ]
        #gc_file_list = [os.path.join(outdir, x) for x in os.listdir(outdir) if x.endswith('gc.ready.gz')]
        gc_hash = module_GCnorm.GC_bias(gc_file_list)
        print(gc_hash)
        ##2.2 GC bias normalization##
        gc_pool = Pool(args.process)
        for gc_file in gc_file_list:
            gc_pool.apply_async(module_GCnorm.GC_norm, args=(gc_file, gc_hash,))
        gc_pool.close()
        gc_pool.join()
        showTime('GC correction done at')
        #os.remove(mpileup_file)
    
    #3# sliding window
    if('3' in args.step):
        cutoff = 0.9
        ##3.1##
        #sliding window
        import module_slidingwindow
        showTime('Sliding window start at')
        gc_norm_list = [ os.path.join(outdir, x+'.gc.norm.gz') for x in out_list ]
        raw_pool = Pool(args.process)
        for norm_file in gc_norm_list:
            raw_pool.apply_async(module_slidingwindow.slidingWindow, 
                                 args=(norm_file, 10000, 2000,cutoff,))
        raw_pool.close()
        raw_pool.join()
        
    #4# filter
    if('4' in args.step):
        cutoff = 0.9
        ##4.1##
        ##get mean and SD ##
        import module_meanFilter
        win_raw_list = [ os.path.join(outdir, x+'.10k.win.raw') for x in out_list ]
        avg, sd = module_meanFilter.getMeanSD(win_raw_list, cent=0.8)
        fhmean = open(outdir+'/mean_sd.report','w')
        print('[Mean]: {}'.format(avg), file= fhmean)
        print('[SD]: {}'.format(sd), file = fhmean)
        
        ##4.2##
        ##filter with mean + 3SD##
        win_pool = Pool(args.process)
        for raw_file in win_raw_list:
            win_pool.apply_async(module_meanFilter.filterwindow, 
                                 args=(raw_file, avg, sd, cutoff,))
        win_pool.close()
        win_pool.join()
    
    #5# smooth
    if('5' in args.step):
        showTime('Smooth start at')
        import module_smooth
        for mychr in out_list:
            win_file = os.path.join(outdir, mychr +'.10k.win')
            md = module_smooth.rtSmooth(win_file, win=10000)
            module_smooth.rpCurve(mychr, md, outdir)


#runRepTime(args.input, args.chr, args.outdir, args.ref)


			
