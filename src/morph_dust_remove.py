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
import sys
sys.path.append(os.path.abspath('/home/janmejoyarch/Dropbox/Janmejoy_SUIT_Dropbox/SUIT_publicity/colormap_project/src'))
from colormap import filterColor

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
    plt.figure(figsize=(10,4), dpi=300)
    plt.subplot(1,2,1)
    plt.imshow(data, origin='lower', vmin=0, cmap= filterColor[ftr_name])
    plt.title("Unfiltered", fontsize=20)
    plt.tick_params(axis='both', labelsize=16)
    plt.subplot(1,2,2)
    plt.imshow(filtered, origin='lower', vmin=0, cmap= filterColor[ftr_name])
    plt.title("Filtered", fontsize=20)
    plt.tick_params(axis='both', labelsize=16)
    plt.savefig(project_path+"reports/dust_removal.pdf")
    plt.show()


def dust_remove(file, thres, plot=None, sav=None):
    hdu= fits.open(file)[0]
    data= hdu.data/flat_data
    header=hdu.header
    h,w= data.shape 
    required_keys= ['CRPIX1', 'CRPIX2', 'R_SUN']
    missing_keys= [key for key in required_keys if key not in header]
    if missing_keys:
        particle_mask= dilation(data<thres, disk(5))
    else:
        col, row, radius= hdu.header['CRPIX1'], hdu.header['CRPIX2'], hdu.header['R_SUN']-50
        mask= np.ones((h,w))*create_circular_mask(h, w, col, row, radius)
        particle_mask= dilation(mask*(data<thres), disk(5))
    filtered=dilation(erosion(data, disk(1)), disk(7))*particle_mask
    dust_grains= data*particle_mask
    filtered_image= (data*np.logical_not(particle_mask))+filtered
    data_crop, filtered_crop= data[y-s:y+s, x-s:x+s], filtered_image[y-s:y+s, x-s:x+s]
    if plot: plotter(data, filtered_image)
    if sav: 
        sav_hdu=fits.PrimaryHDU(filtered_image, header=header)
        sav_hdu.writeto(sav_path, overwrite=True)
        dust_hdu= fits.PrimaryHDU(dust_grains, header=header)
        dust_hdu.writeto(os.path.join(os.path.dirname(sav_path), f'dust_profile_{os.path.basename(sav_path)}'), overwrite=True)
        print(f'File written to \n{sav}')
    

if __name__=='__main__':
    project_path= os.path.expanduser('~/Dropbox/Janmejoy_SUIT_Dropbox/flat_field/morphological_dust_removal_project/')
    filelist= glob.glob(os.path.join(project_path, 'data/raw/*'))
    file= filelist[0]
    ftr_name= file[-9:-5]
    flat_field= '/data1/janmejoy/system_wide_flat_project/products/50_220/BB03_shtr_0_masked_scatter_corrected_reduced_avg_files_flat_lvl1.fits'
    flat_data= fits.open(flat_field)[0].data
    thres=0.97
    sav_path= os.path.join(project_path, 'products', os.path.basename(file))
    x,y,s= 1366, 1393, 100

    dust_remove(file, thres, plot= True, sav=False)
    

    
    
    
    
