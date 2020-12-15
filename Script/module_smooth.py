from matplotlib import pyplot as plt 
import numpy as np
from csaps import csaps
from sklearn import preprocessing

def rtSmooth(file_path, win=10000, norm=True):
    d = np.loadtxt(file_path)
    if norm:
        d[:,1] = preprocessing.scale(d[:,1])
    tab_dp = []
    #
    gap_tol = 20 * win
    min_interval = 100
    roughness = 1e-17
    
    prev_idx = 1
    
    for i in range(1,len(d)):
        prev_pos = d[i-1,0]
        if (d[i,0] > prev_pos + win + gap_tol):
            if i-1-prev_idx > min_interval:
                p = d[prev_idx:i-1,:]
                #f = csaps(p(:,1), p(:,4), roughness, [], p(:,3));
                #cs = CubicSpline(p[:,0], p[:,3])
                cs = csaps(p[:,0], p[:,1], weights=p[:,2] , smooth=roughness)
                
                for j in range(p.shape[0]):
                    pos = int(p[j,0])
                    repTime = round(p[j,1],4)
                    repTime_sm = np.round(cs(pos),4)
                    tab_dp.append([pos, repTime, repTime_sm])
            prev_idx = i
    if i-1-prev_idx > min_interval:
        p = d[prev_idx:i-1,:]
        cs = csaps(p[:,0], p[:,1], weights=p[:,2] , smooth=roughness)
        for j in range(p.shape[0]):
            pos = int(p[j,0])
            repTime = round(p[j,1],4)
            repTime_sm = np.round(cs(pos),4)
            tab_dp.append([pos, repTime, repTime_sm])
    #tab_dp = np.array(tab_dp)
    return(tab_dp)

def rpCurve(mychr, md, outdir):
    md = np.array(md)
    x = md[:,0]/1000000
    y_o = md[:,1]
    y_s = md[:,2]
    plt.figure(figsize=(6,2),dpi=300)
    plt.scatter(x, y_o, s=0.1)
    plt.scatter(x, y_s, s=0.3, color='black')
    plt.xlabel('Chromosomal coordinate (Mb)')
    plt.ylabel('Replication time')
    plt.title(mychr)
    plt.savefig(outdir+'/'+ mychr + '_rpcurve.pdf')