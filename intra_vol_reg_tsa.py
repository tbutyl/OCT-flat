# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 16:35:28 2019

@author: ebm

"""

#import SimpleITK as sitk
from pystackreg import StackReg
import skimage.io as io
#from skimage.filters import gaussian as gauss
from tkinter.filedialog import askdirectory
from sys import exit as ex
from pathlib import Path
import time
#import os
#from skimage.util import img_as_uint
import total_tsa_flat as flat

def sort_key(pth):
    
    #sorts by volume index. Name will be ../stackXX.tif
    #this returns just XX, any length string of numbers thats
    #converted to int
    return int(pth.stem[5:])

def process(pth):
    
    print(f'\nProcessing {pth}')
    
    save_pth = pth / 'intra_reg_stacks'
    
    try:
        save_pth.mkdir()
        
    except FileExistsError:
        ex('Save File for reg stacks already exists. Delete and re-run.')
        
    #Intravolume registration of the fixed volume
    
    print('Doing intravolume registration.')
    
    files = sorted(list(pth.glob('*.tif')), key=sort_key)
    #files = list(pth.glob('*.tif'))
    loop_starts  = time.time()
    
    for file in files:
        start_time = time.time()
        print(f'Processing: {file}')
        fixed = io.imread(str(file))
        fixed = flat.fix(fixed[:,:,:512], divisor=4)
        save_name = save_pth / f'intra_reg_{file.name}'
        io.imsave(arr=fixed.astype('float32'), fname=str(save_name))
        end_time = time.time()
        print(f'{file} was processed in {(end_time - start_time):.2f}s. \
              \n{((end_time - loop_starts)/60):.2f} minutes have elapsed.')
        
def main():
    source = askdirectory()
    if source == '':
            ex("\n\nExited: No file path selected\n\n")
    #sorting(Path(os.path.abspath(source)))
    process(Path(source))

    print('\n\n\n\t\t\t\t---------\n\t\t\t\tCompleted\n\t\t\t\t---------')

if __name__ == '__main__':
    main()