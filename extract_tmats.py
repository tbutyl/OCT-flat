# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 16:42:20 2019

@author: Lab
"""

import numpy as np
from pathlib import Path

pth1 = Path(r'C:\Users\Lab\Desktop\Eric\kv2.1-tsa-week24\10x\2019-12-03_120906\stacks\intra_reg_stacks\vol_transformation_matrices')
#pth2 = Path(r'C:\Users\Lab\Desktop\Eric\kv2.1-tsa-week22\2019-11-19_113709\stacks\intra_reg_stacks\old_attempts\trans_AdvMSq_prev\vol_transformation_matrices')
#pth3 = Path(r'C:\Users\Lab\Desktop\Eric\kv2.1-tsa-week22\2019-11-19_113709\stacks\intra_reg_stacks\old_attempts\transl_gaussreg_first\vol_transformation_matrices')

def sort_key(pth):
    
    #sorts by volume index. Name will be ../tmat_intra_reg_stackXX.tif
    #this returns just XX, any length string of numbers thats
    #converted to int
    return int(pth.stem[20:])

def extract_shifts(pth):
    files = sorted(list(pth.glob('*.txt')), key=sort_key)
    
    mats = []
    for file in files:
        print(file)
        with open(file) as f:
            i=1
            for line in f.readlines():
                if i==3:
                    mats.append(line[21:-4].split())
                    break
                i+=1
    tmats = np.array(mats).astype('float32')

    return tmats

def mag(tmats):
    
    return np.sqrt(tmats*tmats).sum(axis=1)

tmats1 = extract_shifts(pth1)
#tmats2 = extract_shifts(pth2)
#tmats3 = extract_shifts(pth3)
mags1 = mag(tmats1)
#mags2 = mag(tmats2)
#mags3 = mag(tmats3)
