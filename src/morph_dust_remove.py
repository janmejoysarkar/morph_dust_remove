#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 12:47:14 2024
Dust removal with morphology.
2024-06-30: Dilation Radius reduced from 10 px to 6 px. Screenshots taken for report.
2024-07-03: added provision to save the replaced pixels as a separate image.
Changed mask generation kernel size from 6 px to 3 px.
@author: janmejoyarch
"""
import os, glob
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

def plotter(data, filtered):
    fig, ax= plt.subplots(1,2)
    ax[0].imshow(data, origin='lower', vmin=0)
    ax[0].set_title("Unfiltered")
    ax[1].imshow(filtered, origin='lower', vmin=0)
    ax[1].set_title("Filtered")
    plt.show()

def dust_remove(file, thres, plot=None, sav=None):
    hdu= fits.open(file)[0]
    data= hdu.data
    header=hdu.header
    h,w= data.shape 
    required_keys= ['CRPIX1', 'CRPIX2', 'R_SUN']
    missing_keys= [key for key in required_keys if key not in header]
    if missing_keys:
        particle_mask= dilation(data<thres, disk(3))
    else:
        col, row, radius= hdu.header['CRPIX1'], hdu.header['CRPIX2'], hdu.header['R_SUN']-50
        mask= np.ones((h,w))*create_circular_mask(h, w, col, row, radius)
        particle_mask= dilation(mask*(data<thres), disk(3))
    filtered=dilation(erosion(data, disk(1)), disk(6))*particle_mask
    dust_grains= data*particle_mask
    filtered_image= (data*np.logical_not(particle_mask))+filtered
    data_crop, filtered_crop= data[y-s:y+s, x-s:x+s], filtered_image[y-s:y+s, x-s:x+s]
    if plot: plotter(data_crop, filtered_crop)
    if sav: 
        sav_hdu=fits.PrimaryHDU(filtered_image, header=header)
        sav_hdu.writeto(sav_path, overwrite=True)
        dust_hdu= fits.PrimaryHDU(dust_grains, header=header)
        dust_hdu.writeto(os.path.join(os.path.dirname(sav_path), f'dust_profile_{os.path.basename(sav_path)}'), overwrite=True)
        print(f'File written to \n{sav}')
    

if __name__=='__main__':
    project_path= os.path.expanduser('~/Dropbox/Janmejoy_SUIT_Dropbox/flat_field/morphological_dust_removal_project/')
    filelist= glob.glob(os.path.join(project_path, 'data/raw/*'))
    file= filelist[1]
    thres=0.97
    sav_path= os.path.join(project_path, 'products', os.path.basename(file))
    x,y,s= 1366, 1393, 30

    dust_remove(file, thres, plot= True, sav=False)
    

    
    
    
    