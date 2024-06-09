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
    if sav != None: 
        fits.writeto(sav, data_new, overwrite=True)
        print(f'File written to \n{sav}')

if __name__=='__main__':
    file= '/home/janmejoyarch/sftp_drive/suit_data/level1fits/2024/05/14/engg4/SUT_UNP_9999_999999_Lev1.0_2024-05-14T06.35.33.982_4081NB05.fits'
    thres=3500
    sav= '/home/janmejoyarch/Desktop/test.fits'
    dust_remove(file, thres)