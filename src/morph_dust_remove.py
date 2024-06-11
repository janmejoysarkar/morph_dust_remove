#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 12:47:14 2024
Dust removal with morphology.
@author: janmejoyarch
"""
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from skimage.morphology import dilation, disk, erosion

def create_circular_mask(h, w, col, row, radius): 
    '''
    *** creates circular mask of desired size ***
    - h, w: height and width of canvas
    - col, row: column and row of circle center
    - radius= radius of circle
    '''
    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - col)**2 + (Y-row)**2)
    mask = dist_from_center <= radius
    return mask #can be circular mask of any size

def dust_remove(file, thres, sav=None):
    hdu= fits.open(file)[0]
    data= hdu.data
    header=hdu.header
    h,w= data.shape 
    col, row, radius= hdu.header['CRPIX1'], hdu.header['CRPIX2'], hdu.header['R_SUN']-50
    mask= np.ones((h,w))*create_circular_mask(h, w, col, row, radius)
    particle_mask= dilation(mask*(data<thres), disk(10))
    filtered=dilation(erosion(data, disk(1)), disk(6))*particle_mask
    data_new= (data*np.logical_not(particle_mask))+filtered
    fig, ax= plt.subplots(1,2)
    ax[0].imshow(data, origin='lower'); ax[0].set_title("Unfiltered")
    ax[1].imshow(data_new, origin='lower'); ax[1].set_title("Filtered")
    plt.show()
    if (sav == True): 
        sav_hdu=fits.PrimaryHDU(data_new, header=header)
        sav_hdu.writeto(sav_path, overwrite=True)
        print(f'File written to \n{sav}')

if __name__=='__main__':
    file= '/data1/janmejoy/morphological_dust_removal_project/products/SUT_C24_0302_000379_Lev0.5_2024-05-17T05.24.11.637_4081NB08.fits'
    thres=1000
    sav_path= '/data1/janmejoy/morphological_dust_removal_project/products/SUT_C24_0302_000379_Lev0.5_2024-05-17T05.24.11.637_4081NB08.fits'
    dust_remove(file, thres, sav=True)