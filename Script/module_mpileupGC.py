# -*- coding: utf-8 -*-
"""
@author: WANG
"""

###function###
import sys
import gzip
import time
import re
import os

binpath = sys.path[0]

def showTime(text):
    localtime = time.asctime( time.localtime(time.time()) )
    print("{0}: {1}\n".format(text, localtime), file=sys.stderr, flush=True)

def getGcContentForEachBase(chromosome, region_list, wing, myfh):
    window = 2.0 * wing
    if len(region_list) > window:
        G_C = 0
        for i, bs in enumerate(region_list[0:wing]):
            if bs[2].lower() in ('g', 'c'):
                G_C += 1
        for i, bs in enumerate(region_list[wing+1: 2*wing+1]):
            if bs[2].lower() in ('g', 'c'):
                G_C += 1
        for idx, bs in enumerate(region_list[wing: len(region_list) - wing - 1]):
            i = idx + wing
            print_line = [str(x) for x in [chromosome, bs[0], bs[1], bs[2], G_C/window]]
            print("\t".join(print_line),file=myfh)
            offset = 0
            if bs[2].lower() in ('g', 'c'):
                offset += 1
            if region_list[i+1][2].lower() in ('g', 'c'):
                offset -= 1
            if region_list[i - wing][2].lower() in ('g', 'c'):
                offset -= 1
            if region_list[i + wing + 1][2].lower() in ('g', 'c'):
                offset += 1
            G_C += offset
        bs = region_list[len(region_list) - wing - 1]
        print_line = [str(x) for x in [chromosome, bs[0], bs[1], bs[2], G_C/window]]
        print("\t".join(print_line),file=myfh)

#def runmpileupGC(file,wing, uniq_region):
def runmpileupGC(bam, ref, chr, outdir, wing, unique=False,
                 uniq_dir =  binpath+'/hg19.uniquely.mappable.regions',
                 samtools = binpath+'/samtools'):
    chr_id = re.match(r'.*?(\d+)',chr).group(1)
    mychr = 'chr' + chr_id
    #1#
    ##samtools mpileup##
    showTime(mychr+" samtools mpileup started at")

    mpileup_file = "{0}/{1}.mpileup.txt.gz".format(outdir, mychr)
    cmd = "{0} mpileup -f {1} -r {2} {3} |gzip > {4}".format(samtools, ref, chr, bam, mpileup_file)
    print(cmd)
    os.system(cmd)
    
    #2#
    ##2.1 uniquely.mappable.region##
    showTime(mychr+" GC content started at")
    
    regions = []
    if unique:
        uniq_region = uniq_dir+'/'+mychr+'.uniquely.mappable.region'
        fh_uniq = open(uniq_region, 'r')
        for line in fh_uniq:
            sp = line.rstrip().split()
            regions.append([int(x) for x in sp[1:3]])
        fh_uniq.close()
    else:
        regions.append([1,3000000000])
    idx = 0
    st, ed = regions[idx]
    
    ##2.2 GC content##
    gc_file = "{0}/{1}.gc.ready.gz".format(outdir, mychr)
    fho = gzip.open(gc_file, 'wt')
    last_pos = None
    cons_region = []  # consecutive region
    fh_pile = gzip.open(mpileup_file,'rt')
    for line in fh_pile:
        chrom, pos_str, base, dp_str = line.rstrip().split()[:4]
        chr_id = re.match(r'.*?(\d+)',chrom).group(1)
        chrom = 'chr' + chr_id
        pos = int(pos_str)
        
        if st <= pos <= ed:
            #need this pos
            if not last_pos:
                last_pos = pos - 1
            if pos > last_pos + 1:
                getGcContentForEachBase(chrom, cons_region, wing, fho)
                cons_region = []
            cons_region.append([pos, dp_str, base])
            last_pos = pos
        elif pos > ed:
            idx += 1
            try:
                st, ed = regions[idx]
            except IndexError:
                break            
    getGcContentForEachBase(chrom, cons_region, wing, fho)
    fh_pile.close()
    os.remove(mpileup_file)
    showTime(mychr+" GC content done at")

  
