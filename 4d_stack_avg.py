# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 15:02:27 2019

@author: Lab
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 13:37:54 2019

@author: ebm
"""

from tkinter.filedialog import askdirectory
from skimage import io as io
import sys
from pathlib import Path

def sort_key(pth):
    
    return int(pth.stem[9:])

def avg(pth):
    
    save_pth = pth/'4d_avg_stack'
    save_pth.mkdir()
    
    files = sorted(list(pth.glob('*.tif')), key=sort_key)
    #ASSUMING 32 VOLUMES CAPTURED PER FLASH
    divisor = len(files)/32
    
    try:
        #to make sure that the divisor is a multiple of 32
        assert divisor/int(divisor)==1
        divisor = int(divisor)
    except AssertionError:
        sys.exit('The number of stacks was not a multiple of 32.')
        
        
    for i,file in enumerate(files):
        print(file)
        
        stk = io.imread(str(file))
        if i==0:
            avg_stack = np.empty((32,stk.shape[0],stk.shape[1], stk.shape[2]))
        if i < 32:
            avg_stack[i, :, :, :] = stk/divisor
        else:
            avg_stack[i%32,:,:,:]+=stk/divisor
        
    io.imsave(fname=str(save_pth/'avg_timeseries_stack.tif'), arr=avg_stack.astype('float32'))
        



def main():
    source = askdirectory()
    if source == '':
            sys.exit("\n\nExited: No file path selected\n\n")
    #sorting(Path(os.path.abspath(source)))
    avg(Path(source))
    
    print('done')
    

if __name__ == '__main__':
    main()
    