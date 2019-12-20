# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 15:23:07 2019

@author: ebm

processing ORGs
"""
from tkinter.filedialog import askdirectory
import numpy as np
from skimage import io as io
import sys
from pathlib import Path
from pystackreg import StackReg

def sort_key(pth):
    """
    When loading files from the dispersion compensation program, they are named 0-255.tif. When globing to get the list of files to load them, they are ordered not numerically, but as strings.
    This function returns the number in that file name and is used as the key function in the built-in "sorted" function.
    
    This function is important. Passing pth.stem directly to sorted takes the stem from the 
    DIRECTORY path, not the file. 
    
    Now uses Pathlib to deal with path objects, simplified!
    Path.stem is equivalent to:
    file_end = pth.rsplit(os.sep,1)[1]
    file_number = file_end.rstrip('.tif')
    """
    
    return int(pth.stem)

def im_open(pth):

    """Imports a folder containing a stack of images using scikit image collections"""

    try:
        #Double check to make sure pth is a directory.
        assert pth.is_dir()
        #Get list of tif files in pth and sort based on 'stem' which is the index of the B-scan (i.e. y position)
        #files = sorted(pth.glob('*.tif'), key=sort_key) 
        #load the collection
        #raw_stack = io.imread_collection(files)
        raw_stack = io.imread_collection(str(pth/'*.tif'))
        #turn them into a stack that can be manipulated as a multdimensional array.
        stack = io.collection.concatenate_images(raw_stack)
        
        return stack

    except AssertionError:
        #if the assert fails the program stops. This is a bit heavy handed...
        sys.exit("A non-directory object was given to the __open__ function")
        
 
def sorting(top_path):
    
    sr = StackReg(StackReg.AFFINE)
    
    amp_folders = top_path.glob('**/Linear_Amp_FFT2X')
    #vol_ind allows all set of stacks to be saved in one foder for the TSA-style registration.
    #vol_ind=0
    for folder in amp_folders:
        try:
            assert folder.is_dir()
        except AssertionError:
            continue
        print(folder)
        #make a folder named stacks to save the relevant set of stacks in.
        #for example if the stacks are in ../loc1/"timestamp"/oct/log_amp_fft2x
        #make ../loc1/stacks to save in.
        try:
            stack_pth = folder.parents[2] / 'Blink'
            stack_pth.mkdir()
        #however, to not overwrite the folder, skip it if it exists.
        except FileExistsError:
            pass              
        #try except here because disp comp makes ALL folders and linear will be empty most of the time
        #this skips it.
        try:
            stack = im_open(folder)
        except ValueError:
            continue
        else:
            save_path = stack_pth / 'full_realstack.tif'
            avg_save_path = stack_pth / 'AVG_full_realstack.tif'
            aff = sr.register_transform_stack(stack, reference='first', n_frames=1, moving_average=5)
            avg_aff = np.mean(aff, axis=0)
            io.imsave(str(save_path), aff.astype('float32'))
            io.imsave(str(avg_save_path), avg_aff.astype('float32'))
         



def main():
    source = askdirectory()
    if source == '':
        sys.exit("\n\nExited: No file path selected\n\n")
    #sorting(Path(os.path.abspath(source)))
    sorting(Path(source))

    print('\n\n\n\t\t\t\t---------\n\t\t\t\tCompleted\n\t\t\t\t---------')

if __name__ == '__main__':
    main()