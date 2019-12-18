# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 13:37:54 2019

@author: ebm
"""

from tkinter.filedialog import askdirectory
from skimage import io as io
import sys
from pathlib import Path

def avg(pth):
    
    save_pth = pth/'avg_stack'
    save_pth.mkdir()
    
    files = pth.glob('*.tif')
    for i,file in enumerate(files):
        print(file)
        
        stk = io.imread(str(file))
        if i == 0:
            avg_stack = stk
        else:
            avg_stack+=stk
            
    avg_stack/=i
        
    io.imsave(fname=str(save_pth/'all_average.tif'), arr=avg_stack)
        



def main():
    source = askdirectory()
    if source == '':
            sys.exit("\n\nExited: No file path selected\n\n")
    #sorting(Path(os.path.abspath(source)))
    avg(Path(source))
    
    print('done')
    

if __name__ == '__main__':
    main()
    
