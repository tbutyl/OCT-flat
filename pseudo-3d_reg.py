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
from scipy.ndimage import uniform_filter as ufilt


#if enface reg not working, use A-scan to get vascular patterns? do bsub?
def segment_nfl(img):
    """
    Program that does edge detection. boundary fitting, and segmentation to extract the NFL
    for the enface (Y-shift) registration.
    
    Uses a boxcar filter because it's constant time, though gaussian may have slightly better
    outcomes.
    
    Edge detection with prewitt filter
    
    Uses minima detection to pick out NFL/vitreal border. To do this the RPE/OS border must not be in
    the image, becuase it will be the minimum (can use argmax with np.less and a filter of 70 to get RPE,
    INL and NFL). The volume is cut in half immediately to do this, which also saves a lot of time.
    
    After filtering, the averaged data set is subsampled for fitting, greatly speeding fitting.
    
    The segementation is clunky and I assume there is a better way to do it than with the loop,
    but this isn't too time-costly. 
    """   
    x_len, z_len, y_len = img.shape
    #get a divisor so the x and y are sampled at 128
    #(possibly less if they are <128)
    if x_len>128:
        x_div = int(x_len/128)
    else:
        x_div = 1
    if y_len>128:
        y_div = int(y_len/128)
    else:
        y_div = 1
    
    #Do the edge detection
    uimg = ufilt(img[:,int(z_len/2):, :], size=[64,15,64]) #these might change with changing dimensions
    sub_uimg = uimg[::x_div,:,::y_div]
    edges = flat.vol_edges(sub_uimg)
    #for some reason the first and last column minima are at the top of the image,
    #need to cut these out for fitting and make sure to use 126 as x-range
    mins = np.argmin(edges, axis=1)[:,1:-1]
    
    #Set up the fitting
    x = np.arange(0, x_len, x_div)
    y = np.arange(0, y_len, y_div)
    #set up indices for fitting, cut X bc mins were cut 
    #Note: is this where the reversal happens?
    X,Y = np.meshgrid(x[1:-1], y, copy=False)
    X = X.flatten()
    Y = Y.flatten()
    A = np.array([X*0+1, X, Y, X**2, X**2*Y, X**2*Y**2, Y**2, X*Y**2, X*Y]).T
    Z = mins.flatten()
    
    #least squares fit to get values for the terms matrix, A
    #ret contains residuals sum, rank, and, s. usualy for troubleshooting if needed.
    coeffs, *ret = np.linalg.lstsq(A,Z,rcond=None)
    
    #get the actual values with the full x,y dimensions (i.e. inter and extrapolation due to subsample).
    #add back the offset value that was lost with the volume was cut in half, sp the boundary
    #falls at the right place for the full image.
    nfl_boundary = (flat.calc(x_len, y_len, coeffs)+z_len/2).astype('int')
    
    #Set up segmentation
    x = np.arange(x_len)
    y = np.arange(y_len)
    X,Y = np.meshgrid(x,y,copy=False)
    X = X.flatten()
    Y = Y.flatten()
    
    #Nfl thickness in px
    nfl_thickness = 20
    #set up segmented array
    holder = np.empty((x_len, nfl_thickness, y_len))
    #fill array iteratively. May be a more elegant way.
    #if the nfl is close to the image bottom then the fit may include the edge, which cannot be index.
    #this just sets any index that is too large to the edge of the image
    nfl_boundary[np.where(nfl_boundary>z_len-1)[0]] = z_len-1
    #take the nfl layers and put in frames -> i.e. flatten!      
    for z_step in range(nfl_thickness):
        holder[Y, z_step, X] = img[Y,nfl_boundary-z_step, X] #not sure why X and Y are flipped compared to fitting
    #return the averaged NFL for registration
    return np.mean(holder, axis=1)
    
    
def sort_key(pth):
    
    #sorts by volume index. Name will be ../stackXX.tif
    #this returns just XX, any length string of numbers thats
    #converted to int
    return int(pth.stem[5:])


def reg(fix,fix_enface, mov):
      
    #assumes (x,z,y) shape. x and y are same shape
    sr = StackReg(StackReg.TRANSLATION) #may need to update to rigid body
    
    #first do enface by rotating and averaging down stack
    #this will find x-y shift and just the y-shift will be taken
    #rot_fix_mean = np.mean(np.rot90(fix, axes=(0,1)), axis=0)
    mov_enface = segment_nfl(mov)
    rot_mov = np.rot90(mov, axes=(0,1))
    #rot_mov_mean = np.mean(rot_mov, axis=0)
    
    #do the registration
    y_mat = sr.register(fix_enface, mov_enface)
    #copy matrix for all frames in moving stack
    #probably a more elegant way to do this with reeat or tile, but...
    enface_mat = np.zeros((mov.shape[1], 3,3))
    enface_mat[:,0,0] = 1
    enface_mat[:,1,1] = 1
    enface_mat[:,2,2] = 1
    enface_mat[:,1,-1] = y_mat[1,-1]
    
    y_mov = sr.transform_stack(rot_mov, tmats=enface_mat)
    
    #now do x-z registration
    #rotate back to normal orientation
    temp_mov = np.rot90(y_mov, axes=(1,0))
    
    #do the registration
    xz_mat = np.zeros((mov.shape[0],3,3))
    #sr = StackReg(StackReg.AFFINE)
    for i,fix_frame,mov_frame in zip(range(mov.shape[0]),fix,temp_mov):
        #skip frame if empty with np.all(mov_frame==0)?
        tmat = sr.register(fix_frame, mov_frame)
        xz_mat[i] = tmat
    
    #do the transform
    reg_mov = sr.transform_stack(temp_mov, tmats=xz_mat)
    
    return reg_mov

    
def save(stk, save_name):
    
    io.imsave(arr=stk.astype('float32'), fname=str(save_name))
    
def process(pth):
    
    print(f'\nProcessing {pth}')
    
    save_pth = pth / 'pseudo_reg_stacks'
    
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
                fix_enface = segment_nfl(fixed)
            except FileNotFoundError:
                fixed = io.imread(str(file))[:,:,:512]
                fix_enface = segment_nfl(fixed)
                #cuts if flyback wasnt cut during disp comp
                #fixed = flat.fix(fixed[:,:,:512], divisor=2) #fix drift - can't fit w/ breathing artifacts
                save(fixed, save_name)
            avg_stk = fixed/num
        else:
            mov = io.imread(str(file))[:,:,:512]
            reg_mov = reg(fixed,fix_enface, mov)
            save(reg_mov, save_name)
            avg_stk += reg_mov/num
            
        end_time = time.time()
        print(f'{file} was processed in {(end_time - start_time):.2f}s. \
              \n{((end_time - loop_starts)/60):.2f} minutes have elapsed.')
        
    reg_stk = flat.fit_fix(avg_stk)
    avg_pth = pth.parent / 'tsa_stack.tif'
    unreg_pth = pth.parent / 'unreg_tsa_stack.tif'
    save(reg_stk, avg_pth)
    save(avg_stk, unreg_pth)
    
        
def main():
    source = askdirectory()
    if source == '':
            ex("\n\nExited: No file path selected\n\n")
    #sorting(Path(os.path.abspath(source)))
    process(Path(source))
    

    print('\n\n\n\t\t\t\t---------\n\t\t\t\tCompleted\n\t\t\t\t---------')

if __name__ == '__main__':
    main()
