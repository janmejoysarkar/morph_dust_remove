#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 18:06:27 2024
-Created to remove dust spots from SUIT images.
-Uses LED/ PRNU images to indentify dust spots in the images
and corrects them with median filter.

@author: janmejoyarch
"""

import numpy as np
from astropy.io import fits
import os
from astropy.convolution import convolve, Box2DKernel
from skimage.morphology import dilation, disk
from scipy.ndimage import median_filter
from concurrent.futures import ProcessPoolExecutor

def lighten(image_list):
    '''
    Lighten blend
    image_list= List of 2D numpy arrays
    '''
    lighten_blend=np.zeros(np.shape(image_list[0]))
    for image in image_list:
        for i in range(4096):
            for j in range(4096):
                if image[i,j] > lighten_blend[i,j]:
                    lighten_blend[i,j]=image[i,j]
    return(lighten_blend)

def dust_mask(): #make a mask of all dust from PRNU image
    prnu_path= os.path.expanduser('~/Dropbox/Janmejoy_SUIT_Dropbox/flat_field/LED/onboard_PRNU_project/')
    prnu1= fits.open(prnu_path+'data/processed/2024-05-17_prnu_355_ff.fits')[0].data
    prnu2= fits.open(prnu_path+'data/processed/2024-05-17_prnu_355_aa.fits')[0].data
    lighten_blend= lighten([prnu1, prnu2])
    blurred= convolve(lighten_blend, Box2DKernel(4), normalize_kernel=True)
    mask= dilation(blurred<0.99, disk(3))
    print("Mask Generated")
    return (mask)

def filter_image(image_path): #main process to remove dust spots
    hdu= fits.open(image_path)[0]
    print(f'{image_path[-64:]} read')
    data, header= hdu.data, hdu.header
    mask= dust_mask()
    med= median_filter(data, footprint=disk(10))
    filtered= (data*(np.invert(mask)))+(med*mask)
    if save==True: 
        header['COMMENT']=('Dust Corrected', 'Comment')
        sav_hdu= fits.PrimaryHDU(filtered, header=header)
        sav_hdu.writeto(os.path.join(sav, image_path[-64:]), overwrite=True)
    return(filtered)
        
if __name__=='__main__':
    save=True
    sav= os.path.expanduser('~/Dropbox/Janmejoy_SUIT_Dropbox/flat_field/morphological_dust_removal_project/products')
    img_folder='/home/janmejoyarch/Dropbox/Janmejoy_SUIT_Dropbox/flat_field/morphological_dust_removal_project/data/raw/'
    file_list=['SUT_C24_0302_000379_Lev0.5_2024-05-17T05.24.11.637_4081NB08.fits']
    
    path_list= [os.path.join(img_folder, file) for file in file_list]

    with ProcessPoolExecutor() as executor:
        executor.map(filter_image, path_list)
