# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 11:45:00 2019

@author: Eric Miller

Intravolume registration to flatten stacks in cone blink experiments.

Alternative to volume registration due to the small field of view size, leading to a lack of
features, and a difficulty registering the volumes.
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 14:35:20 2019

@author: Eric Miller

This script registers multiple volumes so they can be averaged together.
Choose a folder that has all the stacks saved in the same folder!!

This is for temporal superaveraging or for the cone-blink experiments.

The registration works well if the area imaged is small, such that the retinal curvature approaches zero.
If the field of view is larger, flattening will be required. 

Note: flattening likely could function as the pre intravolume registration step and would be efficacious
for both large and small fields of view. 

The current flattening program based on cross-correlation of A-scans may works, but the flattening may not
be "clean" enough.

Warning: This only works with stacks that have been saved as multipage tiffs!
The expectation is for blink this will be run after volume_separator.py.

For superaveraging, a script will need to be written that puts the mutipage tiff stacks in the same folder.
Happens in matlab, but is part of the processing pipeline.
"""
import SimpleITK as sitk
from pystackreg import StackReg
import skimage.io as io
from skimage.filters import gaussian as gauss
from tkinter.filedialog import askdirectory
from sys import exit as ex
from pathlib import Path
import time
#import os
from skimage.util import img_as_uint, img_as_float32

def process(pth):
    
    print(f'\nProcessing {pth}')
    
    save_pth = pth / 'reg_stacks'
    #tmat_pth = pth / 'transformation_matrices'
    try:
        save_pth.mkdir()
        #tmat_pth.mkdir()
    except FileExistsError:
        ex('Save File for reg stacks or tmats already exists. Delete and re-run.')
        
    
    #tell pystack reg that we will use a translational transformation
    #there shouldn't be intra-volume rotation or shear (there might be for rod blink)
    sr = StackReg(StackReg.TRANSLATION)
    #register to the first slice and output the transfomation matrix without transforming
    #iterate through the files in the stacks folder
    files = pth.glob('*.tif')
    loop_starts  = time.time()
    
    for i,file in enumerate(files):
        #start_time = time.time()
        print(f'Processing: {file}')
        #Intravolume registration of the fixed volume
    
        #load first stack and do a 3d gaussian blur with a sigma=1
        #str needed for imread, otherwise only one frame is loaded
        fixed = io.imread(str(file)) #had gauss(), testing
        
        t_mats = sr.register_stack(gauss(fixed), reference='first', n_frames=1)
        #remove the x shift from all the matrices - horizontal movement isn't relevant here, 
        #the volume should be acquired quickly enough that this isn't a problem.
        t_mats[:,0,-1] = 0
        fixed = sr.transform_stack(fixed, tmats=t_mats)
        
        #Using previous to see if I could get rid of wigge. Didn't seem to work.
        #t_mats = sr.register_stack(gauss(fixed), reference='previous')
        #t_mats[:,1,-1] = 0 
        #fixed = sr.transform_stack(fixed, tmats=t_mats)
        #save the register fixed volume in the parent directory for reference
        save_name = save_pth / f'reg_{file.name}'
        #io.imsave(arr=img_as_uint(fixed), fname=str(save_name))
        io.imsave(arr=img_as_float32(fixed), fname=str(save_name))
        #get fixed out of memory - may not be worth the time?
        print(f'Intravolume registration complete. File saved at {save_name}')
        
        #end_time = time.time()
        #print(f'{file} was processed in {(end_time - start_time):.2f}s. \
              #\n{((end_time - loop_starts)/60):.2f} minutes have elapsed.')
        
        #de;ete emumerate
        #if i==4:
            #ex('auto break')
    end_time = time.time()
    print(f'Run took {(end_time-loop_starts)/60:.2f} minutes')

def main():
    source = askdirectory()
    if source == '':
            ex("\n\nExited: No file path selected\n\n")
    #sorting(Path(os.path.abspath(source)))
    process(Path(source))

    print('\n\n\n\t\t\t\t---------\n\t\t\t\tCompleted\n\t\t\t\t---------')

if __name__ == '__main__':
    main()