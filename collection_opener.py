# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 13:16:06 2019

@author: Lab
"""

from tkinter.filedialog import askdirectory
from skimage import io as io
import sys
from pathlib import Path
from pystackreg import StackReg
from skimage.filters import gaussian as gauss


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
        raw_stack = io.imread_collection(str(pth/'*.tif'))
        #turn them into a stack that can be manipulated as a multdimensional array.
        stack = io.collection.concatenate_images(raw_stack)
        
        return stack

    except AssertionError:
        #if the assert fails the program stops. This is a bit heavy handed...
        sys.exit("A non-directory object was given to the __open__ function")
  

def reg(stack):

    sr = StackReg(StackReg.TRANSLATION)
    t_mats = sr.register_stack(gauss(stack), reference='first', n_frames=1)
    t_mats[:,0,-1] = 0
    reg_stack = sr.transform_stack(stack, tmats=t_mats)  

    return reg_stack    
        
def main():
    source = askdirectory()
    if source == '':
            sys.exit("\n\nExited: No file path selected\n\n")
    #sorting(Path(os.path.abspath(source)))
    stack = im_open(Path(source))
    #reg_stack = reg(stack)
    
    #return reg_stack
    return stack
    


    print('\n\n\n\t\t\t\t---------\n\t\t\t\tCompleted\n\t\t\t\t---------')
"""
if __name__ == '__main__':
    main()
"""