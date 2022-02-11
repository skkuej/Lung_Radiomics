
/*******************************************************
 * Copyright (C) 2019-2022 Eunjin Kim, Hwanho Cho <dmswlskim970606@gmail.com, nara9313@gmail.com>
 * 
 * This file is part of Radiomics Based on Lung Adenocarcinoma for Survival Analysis Project.
 * 
 * Radiomics Based on Lung Adenocarcinoma for Survival Analysis Project can not be copied and/or distributed without the express
 * permission of EJK and HHC
 *******************************************************/

#!/usr/bin/env python
# coding: utf-8
 
## Prepare the library for create peri-tumoral ROIs
import os
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import pandas as pd
import SimpleITK as sitk
from radiomics import imageoperations, glszm, featureextractor, getTestCase
from skimage.morphology import disk, dilation,ball,erosion
from skimage import io

## Load the data lists
pt_path = r'/home/data_path/'
pt_list = os.listdir(pt_path)

imgs,rois = [], []
resol,disk_size = [], []
dilation_3d_arr, erosion_3d_arr, sub_arr = [], [], []

for pt in pt_list:
    ############ Data load and Confirm the Resolution of image ############
    print("\n number %s patient processing..." % pt)
    img=nib.load(pt_path+pt+'\\'+pt+'ct_image.img')
    resol.append(img.header.get_zooms())

    # if x, y resolution are the same,
    if resol[-1][0] == resol[-1][1]: 
        ## We create a disk for dilation and erosion which of radius is 5mm, respectively.
        disk_size.append(np.round(5 / (np.sqrt(resol[-1][0]**2 + resol[-1][1]**2 + resol[-1][2]**2))))
    else : 
        print('%s patient doesn\'t have same x, y resolution %f %f'% (pt, resol[-1][0], resol[-1][1]))
        continue

    img = img.get_fdata()
    imgs.append(img)
    roi = nib.load(pt_path+pt+'\\flip_roi1.nii')
    roi_raw = roi
    roi = roi.get_fdata()
    
    ## Make ROI to be binary values
    roi_min = np.min(roi)
    roi[roi>roi_min] = 1
    roi[roi<=roi_min] = 0
    rois.append(roi)

    
    ############ FOR peri-tumoral ROI, we conduct the process of dilation & erosion and substract ##############
    dilation_3d_arr.append(dilation(roi,selem=ball(disk_size[-1])))
    erosion_3d_arr.append(erosion(roi,selem=ball(disk_size[-1])))
    sub_arr.append(dilation_3d_arr[-1] - erosion_3d_arr[-1])

    ## Check the peri-tumoral ROI
    overlap = img * sub_arr[-1]
    title = ['img', 'roi', 'dilated roi 1cm', 'subtract 1cm','overlap 1cm']
    display=[img,roi,dilation_3d_arr[-1],sub_arr[-1],overlap]
    _,_,nz = np.nonzero(roi)
    plt.figure(figsize=(15,7))
    for i in range(len(display)):
        plt.subplot(2,3,i+1)
        plt.imshow(display[i][:,:,nz[0]],cmap='gray')
        plt.colorbar()
        plt.title(title[i])
    plt.show()

    ni_img = nib.Nifti1Image(sub_arr[-1],roi_raw.affine)
    path_save = r'/home/peri_data/'
    nib.save(ni_img,path_save+'\\'+pt+'_1cm.nii.gz')

        
    
  

