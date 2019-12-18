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
import os
from skimage.util import img_as_uint
import total_tsa_flat as flat

def sort_key(pth):
    
    #sorts by volume index. Name will be ../intra_reg_stackXX.tif
    #this returns just XX, any length string of numbers thats
    #converted to int
    return int(pth.stem[15:])


def process(pth):
    
    print(f'\nProcessing {pth}')
    
    save_pth = pth / 'vol_reg_stacks'
    tmat_pth = pth / 'vol_transformation_matrices'
    try:
        save_pth.mkdir()
        tmat_pth.mkdir()
    except FileExistsError:
        ex('Save File for reg stacks or tmats already exists. Delete and re-run.')
        
    #Intravolume registration of the fixed volume
    
    print('Doing intravolume registration.')
    #load first stack and do a 3d gaussian blur with a sigma=1
    #str needed for imread, otherwise only one frame is loaded
    #fixed = io.imread(str(pth/'stack0.tif')) #had gauss(), testing
    #tell pystack reg that we will use a translational transformation
    #there shouldn't be intra-volume rotation or shear (there might be for rod blink)
    #sr = StackReg(StackReg.TRANSLATION)
    #register to the first slice and output the transfomation matrix without transforming
    #t_mats = sr.register_stack(gauss(fixed), reference='first', n_frames=1)
    #remove the x shift from all the matrices - horizontal movement isn't relevant here, 
    #the volume should be acquired quickly enough that this isn't a problem.
    #t_mats[:,0,-1] = 0
    #fixed = sr.transform_stack(fixed, tmats=t_mats)
    #save the register fixed volume in the parent directory for reference
    #fixed = flat.fix(fixed)
    #save_name = pth.parent / 'fixed.tif'
    #io.imsave(arr=img_as_uint(fixed), fname=str(save_name))
    #io.imsave(arr=fixed.astype('float32'), fname=str(save_name))
    #get fixed out of memory - may not be worth the time?
    #del fixed
    #print(f'Intravolume registration complete. File saved at {save_name}')
    
    #Uncomment to do intravol reg for fixed
    #fixed = io.imread(str(pth/'intra_reg_stack0.tif'))
    #fixed = flat.fix(fixed)
    #save_name = pth.parent / 'fixed.tif'
    #io.imsave(arr=fixed.astype('float32'), fname=str(save_name))
    #get fixed out of memory - may not be worth the time?
    #del fixed
    #print(f'Intravolume registration complete. File saved at {save_name}')
    
    #Volume registration
    
    #switch to tmat folder because i dont know how to set the output folder to save 
    #the transformation mats. elastix.Execute() automatically saves
    #TransformParameters.0.txt in the current directory
    os.chdir(tmat_pth)
    elastixImageFilter = sitk.ElastixImageFilter()
    #reload fixed image with sitk so it is in the right format
    #gfixed = sitk.DiscreteGaussian(sitk.ReadImage(str(pth/'intra_reg_stack0.tif')))
    fixed = sitk.ReadImage(str(pth/'intra_reg_stack0.tif'))
    #set up for multivolume registration
    
    #parameterMapVector = sitk.VectorOfParameterMap()
    #parameterMapVector.append(sitk.GetDefaultParameterMap("affine"))
    #parameterMapVector.append(sitk.GetDefaultParameterMap("bspline"))
    #elastixImageFilter.SetParameterMap(parameterMapVector)
    
    elastixImageFilter.SetFixedImage(fixed)
    parameterMap = sitk.GetDefaultParameterMap('translation')
    elastixImageFilter.SetParameterMap(parameterMap)
    #test 11-8-2019
    #parameterMap['Image Sampler'] = ['Full']
    #parameterMap['Metric'] = ['AdvancedNormalizedCorrelation']
    #parameterMap['Optimiser'] = ['FullSearch']
    
    #iterate through the files in the stacks folder
    files = sorted(list(pth.glob('*.tif')), key=sort_key)
    loop_starts  = time.time()
    for i,file in enumerate(files):
        i=0 #use first stack
        #register to previous stack
        if i > 1:
            fixed = sitk.ReadImage(str(save_name))
            #gfixed = sitk.DiscreteGaussian(fixed)
            elastixImageFilter.SetFixedImage(fixed)
        start_time = time.time()
        print(f'Processing: {file}')
        #load the volume
        #Do affine registration so that volume can shear to do the intravolume registration
        #this code is equivalent to the lines below, but allows the export of the
        #transformation matrices
        #moving = sitk.ReadImage(file)
        #reg_volume = sitk.Elastix(fixed, moving, "affine")
        mv = sitk.ReadImage(str(file))
        #gmv = sitk.DiscreteGaussian(mv)
        elastixImageFilter.SetMovingImage(mv)
        elastixImageFilter.Execute()
        tmat = elastixImageFilter.GetTransformParameterMap()
        
        tif = sitk.TransformixImageFilter()
        tif.SetTransformParameterMap(tmat)
        
        tif.SetMovingImage(mv)
        tif.Execute()
        
        reg_volume = tif.GetResultImage()

        #reg_volume = elastixImageFilter.GetResultImage()
        #can get the transform parameters, but have no idea how to save them
        #maybe just through normal python file manipulation?
        #transformParameterMap = elastixImageFilter.GetTransformParameterMap()
        
        #save volume
        save_name = save_pth / f'reg_{file.name}'
        sitk.WriteImage(reg_volume, str(save_name))
        
        #transformation was saved on execute, now rename.
        #tmat_name = tmat_pth / f'tmat_{file.name}'
        os.rename('TransformParameters.0.txt', f'tmat_{file.stem}.txt')
        
        end_time = time.time()
        print(f'{file} was processed in {(end_time - start_time):.2f}s. \
              \n{((end_time - loop_starts)/60):.2f} minutes have elapsed.')
        
        #de;ete emumerate
        #if i==4:
            #ex('auto break')
    

def main():
    source = askdirectory()
    if source == '':
            ex("\n\nExited: No file path selected\n\n")
    #sorting(Path(os.path.abspath(source)))
    process(Path(source))

    print('\n\n\n\t\t\t\t---------\n\t\t\t\tCompleted\n\t\t\t\t---------')

if __name__ == '__main__':
    main()