
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

import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
from radiomics import featureextractor, getTestCase, imageoperations, glszm



def extraction_GLSZMfeatures(image,mask):
    ## Setting for calculate the radiomics features
    kwargs1 = {'binCount':32,  'verbose': True}
    
    ###### Extract GLSZM features #######
    glszmFeatures = glszm.RadiomicsGLSZM(image, mask, **kwargs1) 
    glszmFeatures.enableAllFeatures()
    glszmFeatures.execute()
    df_glszm = pd.DataFrame(glszmFeatures.featureValues,index = [0])
    
    return df_glszm

    

## Define the data path
pt_path=r'/home/data/'
save_root=r'/home/glszm_feature/'

pt_list=os.listdir(pt_path)
results_pd=pd.DataFrame([])
for i in pt_list:
    ## Load the image and ROI data
    img=sitk.ReadImage(save_root+str(i)+'/'+str(i)+'.img')   
    roi=sitk.ReadImage(save_root+str(i)+'/'+str(i)+'roi.nii') # tumoral ROI, peri_tumoral ROI
    roi_np=sitk.GetArrayFromImage(roi)
    
    ## ROI data to be binary 
    roi_np[roi_np>np.min(roi_np)]=1
    roi_np[roi_np<=np.min(roi_np)]=0
    
    ## Get the 16 GLSZM features 
    results=extraction_GLSZMfeatures(img,sitk.GetImageFromArray(roi_np))
    results.index=[i]
    results_pd=pd.concat([results_pd,results],axis=0)

## Save the GLSZM features to a csv file    
results_pd.to_csv(os.path.join(save_root+'/glszm_features.csv'))

