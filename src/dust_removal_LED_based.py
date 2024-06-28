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
import matplotlib.pyplot as plt
import glob

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

def dust_mask(prnu_path): #make a mask of all dust from PRNU image
    prnu= fits.open(prnu_path)[0].data
    blurred= convolve(prnu, Box2DKernel(4), normalize_kernel=True)
    mask= dilation(blurred<0.99, disk(3))
    print("Mask Generated")
    return (mask)

def filter_image(image_path): #main process to remove dust spots
    hdu= fits.open(image_path)[0]
    print(f'{image_path[-64:]} read')
    data, header= hdu.data, hdu.header
    mask= dust_mask(prnu_path)
    med= median_filter(data, footprint=disk(10))
    filtered= (data*(np.invert(mask)))+(med*mask)
    if save==True: 
        header['COMMENT']=('Dust Corrected', 'Comment')
        sav_hdu= fits.PrimaryHDU(filtered, header=header)
        sav_hdu.writeto(os.path.join(sav, image_path[-64:]), overwrite=True)
    return (data, filtered)

def plot(data, filtered):
    plt.figure()
    plt.subplot(1,2,1)
    plt.imshow(data,vmin=0, vmax=5e4, origin='lower')
    plt.title("Raw Image")
    plt.subplot(1,2,2)
    plt.imshow(filtered, vmin=0, vmax=5e4, origin='lower')
    plt.title("Dust Filtered")
    plt.show()
        
if __name__=='__main__':
    save=False
    project_path= os.path.expanduser('~/Dropbox/Janmejoy_SUIT_Dropbox/flat_field/morphological_dust_removal_project/')
    prnu_path= os.path.join(project_path, 'data/external/2023-11-23_prnu_355_common.fits')
    sav= os.path.join(project_path, 'products')
    file_list=glob.glob(os.path.join(project_path, 'data/raw/*'))
    data, filtered= filter_image(file_list[0])
    plot(data, filtered)

    x,y,s= 1366, 1393, 25
    data_crop, filtered_crop= data[y-s:y+s, x-s:x+s], filtered[y-s:y+s, x-s:x+s]
    print(np.sum(data_crop), np.sum(filtered_crop))
    plot(data_crop, filtered_crop)
    
'''
    with ProcessPoolExecutor() as executor:
        executor.map(filter_image, file_list)
'''