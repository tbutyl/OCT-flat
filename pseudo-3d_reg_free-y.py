# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 17:31:02 2019

@author: ebm
"""

from pystackreg import StackReg
import skimage.io as io
#from skimage.filters import gaussian as gauss
from tkinter.filedialog import askdirectory
from sys import exit as ex
from pathlib import Path
import time
import numpy as np
import total_tsa_flat as flat

#if enface reg not working, use A-scan to get vascular patterns? do bsub?


def sort_key(pth):
    
    #sorts by volume index. Name will be ../stackXX.tif
    #this returns just XX, any length string of numbers thats
    #converted to int
    return int(pth.stem[5:])

def reg_frames(fix, mov, sr=StackReg(StackReg.TRANSLATION)):
    
    tmats = np.zeros((mov.shape[0],3,3))
    for i,fix_frame,mov_frame in zip(range(mov.shape[0]),fix,mov):
        #skip frame if empty with np.all(mov_frame==0)?
        tmat = sr.register(fix_frame, mov_frame)
        tmats[i] = tmat
        
    return tmats

def reg(fix,mov):
      
    #assumes (x,z,y) shape. x and y are same shape
    sr = StackReg(StackReg.TRANSLATION) #may need to update to rigid body
    
    #first do enface by rotating and averaging down stack
    #this will find x-y shift and just the y-shift will be taken
    rot_fix = np.rot90(fix, axes=(0,2))
    rot_mov = np.rot90(mov, axes=(0,2))
    y_shifts = reg_frames(rot_fix, rot_mov)
    #in the matrix the x shift is the y shift and the y shift is the z shift
    y_shifts[:,1,-1] = 0
    
    y_mov = sr.transform_stack(rot_mov, tmats=y_shifts)
    
    #now do x-z registration
    #rotate back to normal orientation
    temp_mov = np.rot90(y_mov, axes=(2,0))
    
    xz_shifts = reg_frames(fix, temp_mov)
    
    #do the transform
    reg_mov = sr.transform_stack(temp_mov, tmats=xz_shifts)
    
    return reg_mov

    
def save(stk, save_name):
    
    io.imsave(arr=stk.astype('float32'), fname=str(save_name))
    
def process(pth):
    
    print(f'\nProcessing {pth}')
    
    save_pth = pth / 'pseudo_reg_stacks_freey'
    
    try:
        save_pth.mkdir()
        
    except FileExistsError:
        #ex('Save File for reg stacks already exists. Delete and re-run.')
        print('DANGER DANGER: {save_path} exists. Continuing and overwriting stacks > 0.')
        
    #Intravolume registration of the fixed volume
    
    print('Doing pseudo 3D registration.')
    
    files = sorted(list(pth.glob('*.tif')), key=sort_key)
    num = len(files)
    loop_starts  = time.time()
    
    for i, file in enumerate(files):
        start_time = time.time()
        print(f'Processing: {file}')
        save_name = save_pth / f'pseudo_reg_{file.name}'
        if i==0:
            try:
                fixed = io.imread(str(save_name))
            except FileNotFoundError:
                fixed = io.imread(str(file))[:,:,:512]
                #cuts if flyback wasnt cut during disp comp
                #fixed = flat.fix(fixed[:,:,:512], divisor=2) #fix drift - can't fit w/ breathing artifacts
                save(fixed, save_name)
            avg_stk = fixed/num
        else:
            mov = io.imread(str(file))[:,:,:512]
            reg_mov = reg(fixed, mov)
            save(reg_mov, save_name)
            avg_stk += reg_mov/num
            
        end_time = time.time()
        print(f'{file} was processed in {(end_time - start_time):.2f}s. \
              \n{((end_time - loop_starts)/60):.2f} minutes have elapsed.')
        
    reg_stk = flat.fix(avg_stk, divisor=4)
    avg_pth = pth.parent / 'tsa_stack_fy.tif'
    save(reg_stk, avg_pth)
    
        
def main():
    source = askdirectory()
    if source == '':
            ex("\n\nExited: No file path selected\n\n")
    #sorting(Path(os.path.abspath(source)))
    process(Path(source))
    

    print('\n\n\n\t\t\t\t---------\n\t\t\t\tCompleted\n\t\t\t\t---------')

if __name__ == '__main__':
    main()
