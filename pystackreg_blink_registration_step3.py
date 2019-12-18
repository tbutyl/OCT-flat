# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 14:41:34 2019

@author: Lab
"""
import skimage.io as io
from tkinter.filedialog import askdirectory
from sys import exit as ex
from pathlib import Path
from skimage.util import img_as_float32
import numpy as np

def im_open(pth):

    """Imports a folder containing a stack of images using scikit image collections"""

    try:
        #Double check to make sure pth is a directory.
        assert pth.is_dir()
        #Get list of tif files in pth and sort based on 'stem' which is the index of the B-scan (i.e. y position)
        #files = sorted(pth.glob('*.tif'), key=sort_key) 
        #load the collection
        raw_stack = io.imread_collection(str(pth/'*.tif'))
        #turn them into a stack that can be manipulated as a multdimensional array.
        stack = io.collection.concatenate_images(raw_stack)
        
        return stack

    except AssertionError:
        #if the assert fails the program stops. This is a bit heavy handed...
        ex("A non-directory object was given to the __open__ function")
  

def timeseries(stack):
    stack_shape = stack.shape
    avg_stack = np.empty((32,stack_shape[1], stack_shape[2]))
    for i in range(32):
        avg_frame = np.average(stack[i::32], axis=0)
        avg_stack[i] = avg_frame
    return avg_stack
    
        
def main():
    source = askdirectory()
    if source == '':
            ex("\n\nExited: No file path selected\n\n")
    #sorting(Path(os.path.abspath(source)))
    src = Path(source)
    stack = im_open(src)
    avg_stack = timeseries(stack)
    save_pth = src / 'timeseries'
    save_pth.mkdir()
    io.imsave(arr=img_as_float32(avg_stack), fname = str(save_pth/'timeseries.tif'))
    
    print('\n\n\n\t\t\t\t---------\n\t\t\t\tCompleted\n\t\t\t\t---------')

if __name__ == '__main__':
    main()