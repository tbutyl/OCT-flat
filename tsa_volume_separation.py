# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 15:28:07 2019

@author: Lab
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 14:24:06 2019

@author: Eric Miller

This software is designed to separate OCT volumes that have been saved into a single stack.
It is needed to separated volumes for the faster blink acquisition, which saves all the data into memory
before saving to the computer.

If there are n volumes, the first n frames are the first frames of the n volumes. 
In effect, this is the "deinterleave" command from ImageJ, but processes a large number of volumes
and saves them.

After separation using this script, the will need to be registered. This can be done with
the temporal super averaging (TSA) software.
"""
from tkinter.filedialog import askdirectory
import numpy as np
from skimage import io as io
import glob
import os
import sys
from pathlib import Path


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
     
    """
    Finds all the folders that have had dispersion compensation done on them.
    
    IMPORTANT: For each 5Hz blink/ORG acquistion, the scan will contain multiple volumes, where
    each frame is separated from the next frame by the number of frames in the volume. 
    To say it another way, Frame 1 of Volume 1 will be the first frame in the stack. Frame 2 of 
    Volume 1 will be 33rd in the stack if there are 32 volumes. Frame 1 of Volume 2 is 2nd in the stack.
    
    The frames (B-scans) need to be separated and concatenated in the correct stack, which then needs to 
    be labeled and saved in a separate multipage tiff. Easy to do in numpy.
    
    However, currently to see the cone response we are acquiring every minute for 10+ minutes (i.e. 10+ sets
    of volumes). To register these volumes, they should all be in the same folder. 
    Therefore BEFORE RUNNING THIS, multiple acquisitions of the same area should be put in a separate folder. 
    This should already have been done for batch dispersion compensation, but that code is ill-commented.
    
    Finally, this can be run on an arbitrarily high set of data, but will do weird things to normal oct,
    so if those are present, you will need to run on the folder that contains each set of repeated acquisitions.
    
    After this code, run the TSA code.
    Note to self: figure out where the TSA code is. ¯\_(ツ)_/¯
    """
    
    amp_folders = top_path.glob('**/Log_Amp_FFT2X')
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
            stack_pth = folder.parents[1] / 'stacks'
            stack_pth.mkdir()
        #however, to not overwrite the folder, skip it if it exists.
        except FileExistsError:
            pass              
        #open
        stack = im_open(folder)
        index = folder.parent.stem.split('_')[0][3:]
        save_path = stack_pth / 'stack{}.tif'.format(int(index)-1)
        io.imsave(str(save_path), stack)
        #vol_ind+=1  

def main():
    source = askdirectory()
    if source == '':
            sys.exit("\n\nExited: No file path selected\n\n")
    #sorting(Path(os.path.abspath(source)))
    sorting(Path(source))

    print('\n\n\n\t\t\t\t---------\n\t\t\t\tCompleted\n\t\t\t\t---------')

if __name__ == '__main__':
    main()