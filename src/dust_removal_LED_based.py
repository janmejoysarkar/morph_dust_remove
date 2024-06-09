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

def dust_mask():
    prnu_path= os.path.expanduser('~/Dropbox/Janmejoy_SUIT_Dropbox/flat_field/LED/onboard_PRNU_project/')
    prnu1= fits.open(prnu_path+'data/processed/2024-02-23_prnu_355_ff.fits')[0].data
    prnu2= fits.open(prnu_path+'data/processed/2024-02-23_prnu_355_aa.fits')[0].data
    lighten_blend= lighten([prnu1, prnu2])
    blurred= convolve(lighten_blend, Box2DKernel(4), normalize_kernel=True)
    mask= dilation(blurred<0.99, disk(3))
    return (mask)

def filter_image(image_path):
    hdu= fits.open(image_path)[0]
    data, header= hdu.data, hdu.header
    mask= dust_mask()
    med= median_filter(data, footprint=disk(10))
    filtered= (data*(np.invert(mask)))+(med*mask)
    sav= os.path.expanduser('~/Dropbox/Janmejoy_SUIT_Dropbox/flat_field/dust_removal_project/products')
    if save==True: 
        header['COMMENT']=('Dust Corrected', 'Comment')
        sav_hdu= fits.PrimaryHDU(filtered, header=header)
        sav_hdu.writeto(os.path.join(sav, image_path[-64:]), overwrite=True)
    print(f'{image_path[-64:]} Completed')
        
if __name__=='__main__':
    save=True
    img_folder='/home/janmejoyarch/sftp_drive/suit_data/level1.1fits/2024/03/28/engg4/'
    file_list=['SUT_C24_0262_000314_Lev1.0_2024-03-28T07.00.12.738_4081NB04.fits', 
                'SUT_C24_0262_000314_Lev1.0_2024-03-28T07.00.52.714_4081NB03.fits',
                'SUT_C24_0262_000314_Lev1.0_2024-03-28T07.01.13.826_4081NB02.fits',
                'SUT_C24_0262_000314_Lev1.0_2024-03-28T07.02.12.234_4081NB05.fits',
                'SUT_C24_0262_000314_Lev1.0_2024-03-28T07.02.53.841_4081NB06.fits',
                'SUT_C24_0262_000314_Lev1.0_2024-03-28T07.03.30.607_4081NB07.fits',
                'SUT_C24_0262_000314_Lev1.0_2024-03-28T07.04.12.810_4081NB08.fits',
                'SUT_C24_0262_000314_Lev1.0_2024-03-28T07.04.54.273_4081BB01.fits',
                'SUT_C24_0262_000314_Lev1.0_2024-03-28T07.05.15.394_4081NB01.fits',
                'SUT_C24_0262_000314_Lev1.0_2024-03-28T07.06.18.321_4081BB02.fits',
                'SUT_C24_0262_000314_Lev1.0_2024-03-28T07.06.37.431_4081BB03.fits']
    
    path_list= [os.path.join(img_folder, file) for file in file_list]
    with ProcessPoolExecutor() as executor:
        executor.map(filter_image, path_list)