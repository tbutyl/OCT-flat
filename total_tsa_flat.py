# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 15:01:38 2019

@author: ebm

For intravolume registration and flattening, all with translation.
"""
import numpy as np
from pystackreg import StackReg
from scipy.ndimage import uniform_filter as ufilt
from skimage.filters import prewitt_h as filt


def vol_edges(stk):
    
    edges = np.empty(stk.shape)
    for i, frame in enumerate(stk):
        edges[i] = filt(frame)
    return edges

def calc(x,y, coeffs):
    
    """
    2-D polynomial function that takes x and y values (integers) and an array of ceofficients to compute
    the NFL boundary.
    """
    
    xarr = np.arange(0,x,1)
    yarr = np.arange(0,y,1)
    X,Y = np.meshgrid(xarr,yarr, copy=False)
    X = X.flatten()
    Y = Y.flatten()
    terms_matrix = np.array([X*0+1, X, Y, X**2, X**2*Y, X**2*Y**2, Y**2, X*Y**2, X*Y])
    return coeffs@terms_matrix

#only works if a power of 2 or divisible by 2 and 4
def fix(stk,divisor=1):
    sr = StackReg(StackReg.TRANSLATION)
    
    #Load Stack
    
    #Step 1 - Remove Wiggles
    print('Wiggles')
    tmat = sr.register_stack(stk, reference='previous')
    #clear y shifts
    #tmat[:,1,-1] = 0
    shifts = tmat[:,0,-1]
    frame_ind = np.arange(len(shifts))
    #fit a line to the shifts to remove drift
    coeff = np.polyfit(frame_ind, shifts,3)
    drift_line = np.poly1d(coeff)
    wiggles = shifts-drift_line(frame_ind)
    tmat[:,0,-1] = wiggles
    
    #Step 1 - fix curvature along y - NOT GOOD IF BREATHING ARTIFACTS PRESENT
    print('Curvature Y')
    #polyfit is too smooth if there are breathing artifacts present - only affect y curvature
    #temp = sr.register_stack(stk[::divisor,:,::divisor], reference='first', n_frames=1, moving_average=5)
    #coeffy = np.polyfit(frame_ind[::divisor], temp[:,1,-1],3)
    #y_line = np.poly1d(coeffy)
    #tmat[:,1,-1] = y_line(frame_ind)
    
    #now just register all frames - bumps should be muted somewhat by using the moving average
    #can gaussian filter if too bumpy
    temp = sr.register_stack(stk, reference='first', n_frames=1, moving_average=5)
    tmat[:,1,-1] = temp[:,1,-1]
    
    
    #use the combined matrix to fix the x-wiggles and the z-y curvature
    print('Transform')
    yless_wiggleless = sr.transform_stack(stk, tmats=tmat)
    
    #Step 2 - fix curvature along x
    print('Curvature X')
    temp_rot = np.rot90(yless_wiggleless, k=-1, axes=(0,2))
    temp_stk = np.empty((temp_rot.shape))
    #temp_stk[:256]=temp_rot[256:]
    #temp_stk[256:]=temp_rot[:256]
    half = int(temp_stk.shape[0]/2)
    temp_stk[:half]=temp_rot[half:]
    temp_stk[half:]=temp_rot[:half]
    
    #ind = int(256/divisor) #normally 256 is halfway thru 512
    ind = int(half/divisor)
    temp = sr.register_stack(temp_stk[::divisor,:,::divisor], reference='first', n_frames=1, moving_average=5)
    final_mat = np.empty((temp.shape))
    final_mat[:ind] = temp[ind:]
    final_mat[ind:] = temp[:ind]
    final_mat[:,0,-1] = 0
    coeffx = np.polyfit(frame_ind[::divisor], final_mat[:,1,-1],3)
    x_line = np.poly1d(coeffx)
    final_mat_whole = np.zeros((frame_ind.shape[0], 3, 3))
    final_mat_whole[:,0,0], final_mat_whole[:,1,1], final_mat_whole[:,2,2] = 1,1,1
    final_mat_whole[:,1,-1] = x_line(frame_ind)
    
    #finally fix the x-z curvature
    print('Transform X')
    fixed = sr.transform_stack(temp_rot, tmats=final_mat_whole)
    final = np.rot90(fixed, k=1, axes=(0,2))
    
    return final

def dewiggle(stk):
    
    print('Wiggles')
    sr = StackReg(StackReg.TRANSLATION)
    tmat = sr.register_stack(stk, reference='previous')
    #clear y shifts
    #tmat[:,1,-1] = 0
    shifts = tmat[:,0,-1]
    frame_ind = np.arange(len(shifts))
    #fit a line to the shifts to remove drift
    coeff = np.polyfit(frame_ind, shifts,3)
    drift_line = np.poly1d(coeff)
    wiggles = shifts-drift_line(frame_ind)
    tmat[:,0,-1] = wiggles
    wiggleless = sr.transform_stack(stk, tmats=tmat)
    
    return wiggleless

def fit_fix(stk):
   
    
    #Step 1 - Remove Wiggles

    #wiggleless = dewiggle(stk)
    wiggleless=stk #####
    #polyfit flatten
    # same idea as nfl segmentation but uses the first half of the volume in depth to detect
    # RPE/OS boundary, then flattens by placing into stack
    print('Flattening')
    medStack = ufilt(wiggleless, (64,15,64))
    x_len, z_len, y_len = np.array(wiggleless.shape).astype('int') #x and y are acutally reversed but im lazy and dont want to chang code
    indx, indy = (int(x_len/128), int(y_len/128))
    #indx,indy=1,1 #####
    uarrt = medStack[::indx, :int(z_len/2), ::indy]
    edges = vol_edges(uarrt)
    mins = np.argmin(edges, axis=1)[:, 1:-1]
    
    x = np.arange(0, x_len, indx)
    y = np.arange(0, y_len, indy)[1:-1]
    
    X,Y = np.meshgrid(x,y,copy=False)
    
    X = X.flatten()
    Y = Y.flatten()
    
    A = np.array([X*0+1, X, Y, X**2, X**2*Y, X**2*Y**2, Y**2, X*Y**2, X*Y]).T
    Z = mins.flatten()
    coeffs, *ret = np.linalg.lstsq(A,Z, rcond=None)
    
    poly_shifts = np.round(calc(x_len, y_len, coeffs).reshape(x_len, y_len)).astype('int')
    maxm = poly_shifts.max()
    
    flattened_stack = np.zeros((z_len+int(maxm-poly_shifts.min()), x_len, y_len))
    counter_array = np.tile(np.arange(z_len, dtype='int')[:, np.newaxis, np.newaxis], (1,x_len,y_len))
    shift_arr = np.tile(maxm-poly_shifts[np.newaxis,:,:], (z_len,1,1))
    
    np.put_along_axis(flattened_stack,shift_arr+counter_array, np.swapaxes(wiggleless,1,0), axis=0)    
    
    final = np.rot90(flattened_stack, k=1, axes=(1,0))
    
    return final