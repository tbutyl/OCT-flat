# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 15:31:03 2019

@author: Lab

"""
import SimpleITK as sitk

pmap = sitk.GetDefaultParameterMap('translation')

fixed = sitk.ReadImage(r'C:\Users\Lab\Desktop\Eric\kv2.1-tsa-week22\2019-11-19_113709\stacks\intra_reg_stacks\intra_reg_stack0.tif')
mv = sitk.ReadImage(r'C:\Users\Lab\Desktop\Eric\kv2.1-tsa-week22\2019-11-19_113709\stacks\intra_reg_stacks\intra_reg_stack30.tif')

gfixed = sitk.DiscreteGaussian(fixed)
gmv = sitk.DiscreteGaussian(mv)

eif = sitk.ElastixImageFilter()
eif.SetFixedImage(gfixed)
eif.SetMovingImage(gmv)

eif.SetParameterMap(pmap)
pmap['WriteResultImage'] = ['False']

eif.Execute()

gresult = eif.GetResultImage()
tmat = eif.GetTransformParameterMap()

sitk.WriteImage(gresult, r'C:\Users\Lab\Desktop\Eric\elastix_g_test\gaussreg_g30.tif')
#sitk.WriteParameterFile(tpm, r'C:\Users\Lab\Desktop\Eric\elastix_g_test\gtmat.txt') #Doesn't work?
#WORKS TO SAVE TMAT
eif.WriteParameterFile(eif.GetTransformParameterMap()[0], r'C:\Users\Lab\Desktop\Eric\elastix_g_test\tmat.txt')

#transformix
tif = sitk.TransformixImageFilter()
tif.SetTransformParameterMap(tmat)

tif.SetMovingImage(mv)
tif.Execute()

result = tif.GetResultImage()

sitk.WriteImage(result, r'C:\Users\Lab\Desktop\Eric\elastix_g_test\gaussreg_30.tif')

eif2 = sitk.ElastixImageFilter()
eif2.SetFixedImage(fixed)
eif2.SetMovingImage(mv)

eif2.SetParameterMap(pmap)

eif2.Execute()

base_result = eif2.GetResultImage()
sitk.WriteImage(base_result, r'C:\Users\Lab\Desktop\Eric\elastix_g_test\base_30.tif')