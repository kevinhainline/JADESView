#!/usr/bin/env python
"""\
JADESView v2
Kevin Hainline

Usage: This file runs EAZY on JADES catalog sources, and plots the SED and chi-square surface
       in a GUI, where the user can also specify where mosaic images are for plotting
       thumbnails. 

	   https://github.com/kevinhainline/JADESView
"""
import os
import copy
import ast
import shutil
import sys
import math
import astropy
import numpy as np
import argparse
import matplotlib.pyplot as plt 
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table

import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import h5py 
import eazy 
import eazy.hdf5 

import webbrowser

from multiprocessing import Pool

from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import tarfile

from astropy.visualization import (ZScaleInterval, ImageNormalize, LinearStretch)

# Set up some colors. 
template_color = '#56B4E9'#'#117733'#'green'
NIRC_photometry_color = '#D55E00'#'#882255'#'red'
HST_photometry_color =  '#F786AA'#'#CC6677'#'lightcoral'
chisq_surface_color = '#E69F00'#'#88CCEE'#'blue'
zspec_color = '#332288'#'orange'

output_plot_folder = ''


def flambda_to_fnu(waveang,flambda):
	return (waveang**2.0) / 2.99792458e+18 * flambda

def FluxtoABMag(flux):
	return (-5.0 / 2.0) * np.log10(flux) - 48.60

def ABMagtoFlux(abmag):
	return np.power(10, (abmag + 48.60)/-2.5)

def RawCatalog_To_Fluxes_Errs_RA_DEC(raw_catalog_file, short_filter_list, blended, convolved, which_aperture_phot, uncertainty_to_use):
	# Start by converting the raw output catalog into an input EAZY py catalog 
	
	# Set up the output file location / name
	raw_catalog_file_tag = ''	
	blended_flag = 'input.'
	convolved_flag = ''
	if (blended == True):
		raw_catalog_file_tag = 'blended'	
		blended_flag = 'blended_input.'
	if (convolved == True):
		convolved_flag = 'conv'
	
	# Here's the input file
	blended_convolved_flag = blended_flag + convolved_flag
	if (convolved == True):
		blended_convolved_flag = blended_convolved_flag + '.'
	output_file_name = 'phot_cat.'+which_aperture_phot+'.'+blended_convolved_flag+'fits'
	
	# And now set up the HDU name from the file. 
	if(which_aperture_phot.startswith('KRON')):
		if (convolved == False):
			hdu_name = 'KRON'
		if (convolved == True):
			hdu_name = 'KRON_CONV'
	elif(which_aperture_phot.startswith('CIRC')):
		if (convolved == False):
			hdu_name = 'CIRC'
		if (convolved == True):
			hdu_name = 'CIRC_CONV'	
	else:
		sys.exit('bad aperture phot!')

	raw_catalog = fits.open(raw_catalog_file, memmap = True)
	raw_ID = raw_catalog[hdu_name].data['ID'].astype(int)
	raw_RA_phot = raw_catalog[hdu_name].data['RA']
	raw_DEC_phot = raw_catalog[hdu_name].data['DEC']
	raw_number_pixels = 1.0*(raw_catalog['SIZE'].data['BBOX_XMAX'] - raw_catalog['SIZE'].data['BBOX_XMIN']+1) * (raw_catalog['SIZE'].data['BBOX_YMAX'] - raw_catalog['SIZE'].data['BBOX_YMIN']+1)
	raw_number_pixels_flag_region = 1.0*(raw_catalog['SIZE'].data['BBOX_XMAX'] - raw_catalog['SIZE'].data['BBOX_XMIN']+21) * (raw_catalog['SIZE'].data['BBOX_YMAX'] - raw_catalog['SIZE'].data['BBOX_YMIN']+21)
	no_pixels_detected = np.where(raw_number_pixels == 0)[0]
	raw_number_pixels[no_pixels_detected] = 1

	number_objects = len(raw_ID)

	bad_pixel_fraction_cut = 0.05

	number_filters = len(short_filter_list)

	raw_flux = np.zeros([number_objects, number_filters])
	raw_flux_errors = np.zeros([number_objects, number_filters])

	for q in range(0, len(short_filter_list)):
	
		try:
			raw_flux[:,q] = raw_catalog[hdu_name].data[short_filter_list[q]+'_'+which_aperture_phot]
			raw_flux_errors[:,q] = raw_catalog[hdu_name].data[short_filter_list[q]+'_'+which_aperture_phot+uncertainty_to_use]
		
			zero_values = np.where((raw_catalog['FLAG'].data[short_filter_list[q]+'_FLAG'] == -1) | (raw_catalog['FLAG'].data[short_filter_list[q]+'_FLAG']*1.0/raw_number_pixels_flag_region > bad_pixel_fraction_cut) | (raw_flux_errors[:,q] == 0.00) | np.isnan(raw_flux[:,q]) | np.isnan(raw_flux_errors[:,q]) | np.isinf(raw_flux[:,q]) | np.isinf(raw_flux_errors[:,q]))[0]
	
			raw_flux[zero_values,q] = -9999.
			raw_flux_errors[zero_values,q] = -9999.
		except KeyError:
			raw_flux[:,q] = -9999
			raw_flux_errors[:,q] = -9999

	return raw_ID.astype('int'), raw_flux, raw_flux_errors, raw_RA_phot, raw_DEC_phot

# For an entry redshift, tempfilt, and chi2grid, this returns the closest chisq value. 
def get_chisq(tempfilt_zgrid, chi2fit, z_a_value):
	try:
		return np.round(chi2fit[np.where(tempfilt_zgrid == np.round(z_a_value,2))[0][0]],2)
	except IndexError:
		min_difference = np.min(abs(tempfilt_zgrid - np.round(z_a_value,2)))
		if (min_difference < eazy_self.param['Z_STEP']):
			return np.round(chi2fit[np.argmin(abs(tempfilt_zgrid - np.round(z_a_value,2)))],2)
		else:
			sys.exit("New chisq not found?")

# Given an RA, DEC, and a namestub ("JADES-GS"), this returns the full name of the source
def RADEC_to_RADECName(RA, DEC, namestub):
	RA_string = str(round(RA,5))
	DEC_string =  str(round(DEC,5))
	if (RA > 0):
		RA_string = '+'+RA_string
	if (DEC > 0):
		DEC_string = '+'+DEC_string
	
	ID_RA_DEC =  namestub+RA_string+DEC_string

	return ID_RA_DEC

# This is the function that allows you to get the image cutout given an RA/DEC.
def return_image_cutout(objRA, objDEC, image_hdu, image_wcs, thumbnail_size):

	rasize_value = thumbnail_size#2.0
	decsize_value = thumbnail_size#2.0
	cosdec_center = math.cos(objDEC * 3.141593 / 180.0)

	# Set the position of the object
	position = SkyCoord(str(objRA)+'d '+str(objDEC)+'d', frame='fk5')
	size = u.Quantity((decsize_value, rasize_value), u.arcsec)
	
	image = image_hdu.data
	try:
		image_cutout = Cutout2D(image, position, size, wcs=image_wcs)
	except astropy.nddata.utils.NoOverlapError:
		image_cutout = -9999
	
	return image_cutout#.data	


# This returns the EAZY output values and arrays for plotting
def return_eazy_output(object_ID, eazy_self, zout, primary_z = -9999, alt_z = 0.0):#, h5file = False, savefile = None, savesed = False, savechisq = False):
	
	# First, the redshift grid for the chisq surface. 
	output_zgrid = eazy_self.zgrid#np.arange(self.param['Z_MIN'], self.param['Z_MAX'], self.param['Z_STEP'])

	# the index of the source. 
	i = np.where(eazy_self.OBJID == object_ID)[0][0]

	# Let's get the data from the eazy-py file 
	tempfilt_zgrid = eazy_self.zgrid
	chi2fit = eazy_self.chi2_fit[i,:]

	output_chisq = eazy_self.chi2_fit[i,:]
	
	z_a_value = np.round(zout['z_raw_chi2'][i],2)
	
	if (primary_z < 0):
		data = eazy_self.show_fit(object_ID, show_fnu = True, zshow = z_a_value, get_spec = True)
		output_redshift = z_a_value
	else:
		data = eazy_self.show_fit(object_ID, show_fnu = True, zshow = primary_z, get_spec = True)
		output_redshift = primary_z

	output_wavelength = data['templz']
	output_flux = data['templf']*1e3
	output_phot_wavelength = data['pivot']
	output_phot_template = data['model']*1e3
	output_phot_err_template = data['emodel']*1e3
	output_phot = data['fobs']*1e3
	output_phot_err = data['efobs']*1e3
	output_tef = data['tef']

	return output_wavelength, output_flux, output_phot_wavelength, output_phot_template, output_phot_err_template, output_phot, output_phot_err, output_redshift, tempfilt_zgrid, chi2fit

def return_version_string(input_file_name):
	if (input_file_name.split('/')[-1].startswith('hlsp')):
		try:
			version_string = input_file_name.split('/')[-1].split('_')[-2]
		except IndexError:
			version_string = 'V???'
		
	elif (input_file_name.split('/')[-1].startswith('phot')):
		try:
			version_string = input_file_name.split('/')[-1].split('.')[1]+'.'+input_file_name.split('/')[-1].split('.')[2]+'.'+input_file_name.split('/')[-1].split('.')[3]
		except IndexError:
			version_string = 'V???'
	else:
		version_string = 'V???'

	return version_string

# Here's the Plot GUI information
class PlotGUI:
	def __init__(self, root):
		self.root = root
		self.root.title("JADESView v2")
  		
		self.current_ID = -9999
		self.current_thumbnail_size = -9999
        
		# Main container frame
		main_frame = tk.Frame(root)
		main_frame.pack(fill=tk.BOTH, expand=True)
        
		# Top frame for graph
		plot_frame = tk.Frame(main_frame)
		plot_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
		# Create figure and axis for plot
		self.fig, self.ax = plt.subplots(1,2, figsize=(14, 5))
		self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
		self.canvas_widget = self.canvas.get_tk_widget()
		self.canvas_widget.pack(fill=tk.BOTH, expand=True)
 
		# For clicking.
		self.canvas.mpl_connect("button_press_event", self.on_plot_click)       
		
		if (args_plot_thumbnails):
			# Bottom frame for image grid
			image_grid_frame = tk.Frame(main_frame)
			image_grid_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
			
			# Depending on how many filters you provide, it will change the number of
			# rows, and the number of columns
			if (number_image_filters <= 16):
				n_image_rows = 2
				n_columns = 8
				total_possible = 16
			elif (number_image_filters <= 24):
				n_image_rows = 3
				n_columns = 8
				total_possible = 24
			else:
				n_image_rows = 3
				n_columns = 10
				total_possible = 30
			# This creates the plot that will be appended below the SED.
			self.image_fig, self.image_axes = plt.subplots(n_image_rows, n_columns, figsize=(14, 4))
	
			self.image_fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
			self.image_fig.subplots_adjust(wspace = 0.05, hspace = 0.01)
			self.image_canvas = FigureCanvasTkAgg(self.image_fig, master=image_grid_frame)
			self.image_canvas_widget = self.image_canvas.get_tk_widget()
			self.image_canvas_widget.pack(fill=tk.BOTH, expand=True)
			
			# Load initial images into grid
			# Initialize with zeroed images. 
			self.images = [np.zeros([50,50]) for _ in range(total_possible)]
			for ax, img in zip(self.image_axes.flatten(), self.images):
				ax.imshow(img, cmap='gray', origin = 'lower')
				ax.axis('off')
			self.image_canvas.draw()

        # # # # # #
		# BUTTONS #
		# # # # # #
		 
		# Bottom row of centered buttons
		second_row_frame = tk.Frame(root)
		second_row_frame.pack(side=tk.BOTTOM, fill=tk.X)
 
  		# Update the plot after entering values
		update_button = tk.Button(second_row_frame, text="Update (Enter)", command=self.update_plot, font=("Helvetica", 12, "bold"))
		update_button.pack(side=tk.LEFT)
		#update_button.config(width = 40 )
       	
		# Button to reset the SED plot axes
		reset_button = tk.Button(second_row_frame, text="Reset Axes", command=self.reset_plot, font=("Helvetica", 12, "bold"))
		reset_button.pack(side=tk.LEFT)

		# Button to quit the program
		quit_button = tk.Button(second_row_frame, text="Quit", command=root.quit, font=("Helvetica", 12, "bold"))
		quit_button.pack(side=tk.RIGHT)
 
 		# Button to save the output plt as a png. 
		save_button = tk.Button(second_row_frame, text="Save as PNG", command=self.save_gui, font=("Helvetica", 12, "bold"))
		save_button.pack(side=tk.RIGHT)

 		# Button to save the output plt as a png. 
		save_button = tk.Button(second_row_frame, text="Save EAZY SED + Chisq", command=self.save_eazy_sed, font=("Helvetica", 12, "bold"))
		save_button.pack(side=tk.RIGHT)

		if (args_plot_thumbnails):
	 		# Button to save the output plt as a png. 
			save_button = tk.Button(second_row_frame, text="Save Thumbnails", command=self.save_thumbnails, font=("Helvetica", 12, "bold"))
			save_button.pack(side=tk.RIGHT)

 		# Button to save the output plt as a png. 
		fitsmap_link_button = tk.Button(second_row_frame, text="Fitsmap", command=self.open_fitsmap, font=("Helvetica", 12, "bold"))
		fitsmap_link_button.pack(side=tk.RIGHT)

       
		# Controls frame
		controls = tk.Frame(root)
		controls.pack(side=tk.BOTTOM, fill=tk.X)

		# ID text box
		tk.Label(controls, text="ID").pack(side=tk.LEFT)
		self.id_entry = tk.Entry(controls, width=8)
		self.id_entry.pack(side=tk.LEFT)
		self.id_entry.insert(0, "")
                
        # SED plot Xmin
		tk.Label(controls, text="X Min:").pack(side=tk.LEFT)
		self.xmin_entry = tk.Entry(controls, width=5)
		self.xmin_entry.pack(side=tk.LEFT)
		self.xmin_entry.insert(0, "0")
        
        # SED plot Xmax
		tk.Label(controls, text="X Max:").pack(side=tk.LEFT)
		self.xmax_entry = tk.Entry(controls, width=8)
		self.xmax_entry.pack(side=tk.LEFT)
		self.xmax_entry.insert(0, "5.5")

        # SED plot Ymin
		tk.Label(controls, text="Y Min:").pack(side=tk.LEFT)
		self.ymin_entry = tk.Entry(controls, width=8)
		self.ymin_entry.pack(side=tk.LEFT)
		self.ymin_entry.insert(0, "")
        
        # SED plot Ymax
		tk.Label(controls, text="Y Max:").pack(side=tk.LEFT)
		self.ymax_entry = tk.Entry(controls, width=5)
		self.ymax_entry.pack(side=tk.LEFT)
		self.ymax_entry.insert(0, "")

        # SED plot alternate redshift
		tk.Label(controls, text="Alt z:").pack(side=tk.LEFT)
		self.altz_entry = tk.Entry(controls, width=5)
		self.altz_entry.pack(side=tk.LEFT)
		self.altz_entry.insert(0, "")

  		# Update the plot after entering values
		previous_button = tk.Button(controls, text="Previous", command=self.previous_object)
		previous_button.pack(side=tk.LEFT)
		#update_button.config(width = 40 )

  		# Update the plot after entering values
		next_button = tk.Button(controls, text="Next", command=self.next_object)
		next_button.pack(side=tk.LEFT)
		#update_button.config(width = 40 )

     
    	# Enter an updated chisq ymax value
		self.chisq_ymax_entry = tk.Entry(controls, width=5)
		self.chisq_ymax_entry.pack(side=tk.RIGHT)
		tk.Label(controls, text="Chisq Y Max:").pack(side=tk.RIGHT)
		self.chisq_ymax_entry.insert(0, "")

  		# A button to change the chisq plot to P(z) instead. 
		p_of_z_button = tk.Button(controls, text="P(z)", command=self.p_of_z_plot)
		p_of_z_button.pack(side=tk.RIGHT)

		if (args_plot_thumbnails):
			# Enter an updated chisq ymax value
			self.thumbnail_size = tk.Entry(controls, width=3)
			self.thumbnail_size.pack(side=tk.RIGHT)
			tk.Label(controls, text="Thumbnail Size (arcsec):").pack(side=tk.RIGHT)
			self.thumbnail_size.insert(0, "2")

		# Bind the return button to specify updating the plot
		self.root.bind('<Return>', lambda event: self.update_plot())
		self.root.bind('<KP_Enter>', lambda event: self.update_plot())
		self.root.bind('<Next>', lambda event: self.next_object())
		self.root.bind('<Prior>', lambda event: self.previous_object())

		# Plot the first object in the list. 
		self.id_entry.delete(0, tk.END)  # Delete current text
		self.id_entry.insert(0,str(first_object))
		#self.fig.tight_layout()
		self.update_plot()

    # The main driver function that updates plots in the GUI. 
	def update_plot(self):

		self.ax[0].clear()
		self.ax[1].clear()
		
		# Get the ID entry.
		if (self.id_entry.get() == ''):
			print("NOTE: No ID specified, Returning to the first ID number.")
			self.id_entry.delete(0, tk.END)  # Delete current text
			self.id_entry.insert(0,str(first_object))
		
		object_ID = int(self.id_entry.get())
		
		# Check to make sure the user isn't trying to generate a thumbnail with a negative
		# size.	
		if (args_plot_thumbnails):
			if (self.thumbnail_size.get() == ''):
				print("NOTE: Thumbnail size must be above 0. Resetting to 2 arcseconds.")
				self.thumbnail_size.delete(0, tk.END)
				self.thumbnail_size.insert(0,"2") 
				current_thumbnail_size = 2
			else:
				current_thumbnail_size = float(self.thumbnail_size.get())

				if (current_thumbnail_size <=0):
					print("NOTE: Thumbnail size must be above 0. Resetting to 2 arcseconds.")
					self.thumbnail_size.delete(0, tk.END)
					self.thumbnail_size.insert(0,"2") 
					current_thumbnail_size = 2
		
		# Find the index in the input flux file.
		object_ID_index = np.where(ID_values == object_ID)[0]

		# Make sure that this object ID is actually in the list, otherwise
		# just revert to the current object ID.
		input_eazy_index = np.where(eazy_self.OBJID == object_ID)[0]
		if (len(input_eazy_index) <= 0):
			print("NOTE: ID "+str(object_ID)+" is not in the input list.")
			self.id_entry.delete(0, tk.END)
			self.id_entry.insert(0,str(self.current_ID)) 
			object_ID = int(self.id_entry.get())
			object_ID_index = np.where(ID_values == object_ID)[0][0]
		else:
			object_ID_index = object_ID_index[0]
		
		# Set the RA/Dec for the source.
		objRA = RA_values[object_ID_index]
		objDEC = DEC_values[object_ID_index]

		# And here it's hard coded to specify that the source is
		# in either GOODS-N or GOODS-S
		if ((objDEC > 62) & (objDEC < 63)):
			namestub = 'JADES-GN'
		elif ((objDEC > -28.1) & (objDEC < -27.5)):
			namestub = 'JADES-GS'
		else:
			namestub = 'UNKNOWN'
			
		# Get the JADES RA/DEC ID value	
		JADES_ID = RADEC_to_RADECName(objRA, objDEC, namestub) 
		
		# Here are the colors that are used in the plotting.
		template_color = '#56B4E9'#'#117733'#'green'
		NIRC_photometry_color = '#D55E00'#'#882255'#'red'
		HST_photometry_color =  '#F786AA'#'#CC6677'#'lightcoral'
		chisq_surface_color = '#E69F00'#'#88CCEE'#'blue'
		zspec_color = '#332288'#'orange'
		alternate_color = 'grey'

		# This is the best fit values from EAZY, which we'll always want to keep
		# on the plot.
		bf_output_wavelength, bf_output_flux, bf_output_phot_wavelength, bf_output_phot_template, bf_output_phot_err_template, bf_output_phot, bf_output_phot_err, bf_z_a_value, bf_tempfilt_zgrid, bf_chi2fit = return_eazy_output(object_ID, eazy_self, zout, primary_z = -9999, alt_z = 0.0)

		# This determines whether or not there is an alternate redshift specified
		# by the user. 
		if (self.altz_entry.get() == ''):
			output_wavelength, output_flux, output_phot_wavelength, output_phot_template, output_phot_err_template, output_phot, output_phot_err, z_a_value, tempfilt_zgrid, chi2fit = return_eazy_output(object_ID, eazy_self, zout, primary_z = -9999, alt_z = 0.0)
		else:

			# Make sure that the alternate redshift is between the minimum and maximum values.
			# Otherwise, reset to the best-fit redshift.
			if ((float(self.altz_entry.get()) > eazy_self.param['Z_MIN']) and (float(self.altz_entry.get()) < eazy_self.param['Z_MAX'])):
				output_wavelength, output_flux, output_phot_wavelength, output_phot_template, output_phot_err_template, output_phot, output_phot_err, z_a_value, tempfilt_zgrid, chi2fit = return_eazy_output(object_ID, eazy_self, zout, primary_z = float(self.altz_entry.get()), alt_z = 0.0)
			else:
				self.altz_entry.delete(0, tk.END)
				self.altz_entry.insert(0,"") 
				print("NOTE: Alternate redshift must be between "+str(eazy_self.param['Z_MIN'])+" and "+str(eazy_self.param['Z_MAX']))
				output_wavelength, output_flux, output_phot_wavelength, output_phot_template, output_phot_err_template, output_phot, output_phot_err, z_a_value, tempfilt_zgrid, chi2fit = return_eazy_output(object_ID, eazy_self, zout, primary_z = -9999, alt_z = 0.0)

		# Let's plot the SED:
		pos_flux_errors = np.where(output_phot_err > 0)[0]
		if (self.altz_entry.get() == ''):
			self.ax[0].plot(output_wavelength/1e4, output_flux, color = template_color, lw = 2, zorder = 0, label = '$z_a$ = '+str(round(z_a_value,2)))
		else:
			self.ax[0].plot(bf_output_wavelength/1e4, bf_output_flux, color = 'grey', lw = 2, zorder = 0, label = '$z_a$ = '+str(round(bf_z_a_value,2)), alpha = 0.5)
			self.ax[0].plot(output_wavelength/1e4, output_flux, color = template_color, lw = 2, zorder = 0, label = '$z_{\mathrm{alt}}$ = '+str(round(z_a_value,2)))
		self.ax[0].scatter(output_phot_wavelength/1e4, output_phot_template, marker = 's', s = 110, edgecolor = template_color, color = 'None', zorder = 5)
		
		# And now let's plot the photometry, along with the flux uncertainties, and the filter widths. 
		self.ax[0].scatter(output_phot_wavelength[pos_flux_errors]/1e4, output_phot[pos_flux_errors], color = color_filters[pos_flux_errors], s = 50, zorder = 10)
		self.ax[0].errorbar(output_phot_wavelength[pos_flux_errors]/1e4, output_phot[pos_flux_errors],  xerr = filter_bw[pos_flux_errors]/2.0, yerr = output_phot_err[pos_flux_errors], color = 'black', ls = 'None', alpha = 0.4)

		# Currently, let's look at flux on a logarithmic axis, in the future I'll make
		# this something you can toggle.	
		self.ax[0].semilogy()
		self.ax[0].set_xlabel('Observed Wavelength ($\mu$m)')
		self.ax[0].set_ylabel('Flux (nJy)')
		
		# Make sure that we specify that the SED fit is from EAZY-py.
		if (args_asada_cgm == True):
			self.ax[0].text(0.8, 0.06, 'EAZY-py fit, Asada+24 CGM', fontsize = 12, horizontalalignment='center',
				verticalalignment='center', transform=self.ax[0].transAxes)
		else:
			self.ax[0].text(0.85, 0.06, 'EAZY-py fit', fontsize = 12, horizontalalignment='center',
				verticalalignment='center', transform=self.ax[0].transAxes)
					
		# And now let's figure out the X and Y limits for the plot. 
		# Get the minimum and maximum x value.
		xmin = float(self.xmin_entry.get())
		xmax = float(self.xmax_entry.get())

		self.ax[0].set_xlim(xmin, xmax)

		# Which are the positive fluxes? 
		pos_fluxes = np.where(output_phot > 0)[0]

		if (self.ymin_entry.get() == ''):
			ymin = np.min(output_phot[pos_fluxes]) - (0.4 * np.min(output_phot[pos_fluxes]))
			if (ymin < 0.01):
				ymin = 0.01
		else:
			ymin = float(self.ymin_entry.get())
		if (self.ymax_entry.get() == ''):
			ymax = np.max(output_phot[pos_fluxes]) + (10.0 * np.max(output_phot[pos_fluxes]))
			if (ymax > (100.*np.median(output_phot[pos_fluxes]))):
				ymax = 100.*np.median(output_phot[pos_fluxes])
		else:
			ymax = float(self.ymax_entry.get())

		self.ax[0].set_ylim(ymin, ymax)

		# Here's the title for the plot
		self.ax[0].set_title(JADES_ID+' --- ID '+str(int(object_ID)), fontsize = 15)
	
		# And let's set the legend. 
		self.ax[0].legend(loc = 2, fontsize = 12, frameon=False)

		# First, make sure the user specifies that they want thumbnails
		if (args_plot_thumbnails):
			# We only update the thumbnails if the ID has changed. It takes too long
			# otherwise. 
			if ((object_ID != self.current_ID) or (current_thumbnail_size != self.current_thumbnail_size)):
				self.images = []
				for q in range(0, number_images):
					filter_image_cutout = return_image_cutout(objRA, objDEC, image_hdu_all[q], image_wcs_all[q], float(self.thumbnail_size.get()))
					if (filter_image_cutout == -9999):
						print("Note: No overlap with filter "+str(q))
						self.images.append(np.zeros([50,50]))
					else:
						self.images.append(filter_image_cutout.data)

					image_iterator = np.arange(0, number_images, 1)
					for j, ax, img in zip(image_iterator, self.image_axes.flatten(), self.images):
		
						ax.clear()
						indexerror = 0
						try:
							norm = ImageNormalize(img, interval=ZScaleInterval(), stretch=LinearStretch())
						except IndexError:
							indexerror = 1
						except UnboundLocalError:
							indexerror = 1
										
						if (indexerror == 0):
							ax.imshow(img, cmap='gray', origin = 'lower', aspect='equal', norm = norm)
						else:
							ax.imshow(img, cmap='gray', origin = 'lower', aspect = 'equal')
		
						ax.text(0.51, 0.96, all_images_filter_name[j].split('_')[1], transform=ax.transAxes, fontsize=12, fontweight='bold', ha='center', va='top', color = 'black')
						ax.text(0.5, 0.95, all_images_filter_name[j].split('_')[1], transform=ax.transAxes, fontsize=12, fontweight='bold', ha='center', va='top', color = 'white')
		
						if ((round(SNR_values[object_ID_index,j],2) > -100) and (round(SNR_values[object_ID_index,j],2) < 100)):
							ax.text(0.5, 0.06, 'SNR = '+str(round(SNR_values[object_ID_index,j],2)), transform=ax.transAxes, fontsize=12, fontweight='bold', horizontalalignment='center', color = 'black')
							ax.text(0.5, 0.05, 'SNR = '+str(round(SNR_values[object_ID_index,j],2)), transform=ax.transAxes, fontsize=12, fontweight='bold', horizontalalignment='center', color = 'white')
						elif(round(SNR_values[object_ID_index,j],2) > 100):
							ax.text(0.5, 0.06, 'SNR > 100', transform=ax.transAxes, fontsize=10, fontweight='bold', horizontalalignment='center', color = 'black')
							ax.text(0.5, 0.05, 'SNR > 100', transform=ax.transAxes, fontsize=10, fontweight='bold', horizontalalignment='center', color = 'white')
		
						ax.axis('off')
					self.image_canvas.draw()

		# Here's the current ID. 
		self.current_ID = object_ID
		if (args_plot_thumbnails):
			self.current_thumbnail_size = current_thumbnail_size
		
		# Now, let's plot the chi-square surface. 
		self.ax[1].plot(tempfilt_zgrid, chi2fit, color = chisq_surface_color, zorder = 40)

		# Let's grab the chisq value at the specified redshift. 
		bf_chisq_value = get_chisq(bf_tempfilt_zgrid, bf_chi2fit, bf_z_a_value)
		chisq_value = get_chisq(tempfilt_zgrid, chi2fit, z_a_value)

		# Got to get the index within the EAZY object
		objid_index = np.where(eazy_self.OBJID == int(self.current_ID))[0][0]
		
		# Here are the output redshift values from the EAZY fit
		z160 = zout['z160'][objid_index]
		z840 = zout['z840'][objid_index]
		
		z500 = zout['z500'][objid_index]

		upper_error = z840 - z500
		lower_error = z500 - z160

		z025 = zout['z025'][objid_index]
		z975 = zout['z975'][objid_index]
		
		self.ax[1].set_xlabel(r'z$_{phot}$')
		self.ax[1].set_ylabel(r'$\chi^2$')
		if (self.altz_entry.get() == ''):
			self.ax[1].axvspan(bf_z_a_value-0.05, bf_z_a_value+0.05, color = template_color, label = 'z = '+str(bf_z_a_value)+', $\chi^2 = $'+str(bf_chisq_value), zorder = 20)
		else:
			self.ax[1].axvspan(bf_z_a_value-0.05, bf_z_a_value+0.05, color = 'grey', label = 'z = '+str(bf_z_a_value)+', $\chi^2 = $'+str(bf_chisq_value), zorder = 19)
			self.ax[1].axvspan(z_a_value-0.05, z_a_value+0.05, color = template_color, label = 'z$_{\mathrm{alt}}$ = '+str(z_a_value)+', $\chi^2 = $'+str(chisq_value), zorder = 20)
	
		self.ax[1].axvspan(z500-0.05, z500+0.05, color = 'red', label = 'z500 = '+str(round(z500,2))+'$^{+'+str(round(upper_error,2))+'}$'+'$_{-'+str(round(lower_error,2))+'}$', alpha = 0.2, zorder = 3)

		self.ax[1].axvspan(z025, z975, color = 'grey', zorder = 0, alpha = 0.1)
		self.ax[1].axvspan(z160, z840, color = 'grey', zorder = 1, alpha = 0.05)

		version_string = return_version_string(args_input_file)

		if (default_convolved == True):
			self.ax[1].set_title('Photometry: '+version_string+', '+str(args_aperture)+', Convolved', fontsize = 15)
		else: 
			self.ax[1].set_title('Photometry: '+version_string+', '+str(args_aperture)+', Not Convolved', fontsize = 15)

		
		chisq_legend = self.ax[1].legend(loc = 2, fontsize = 12, frameon=False)
		chisq_legend.set_zorder(50)
		self.ax[1].set_ylim(-10, np.max(chi2fit)+0.1 * np.max(chi2fit))
		if (self.chisq_ymax_entry.get() != ''):
			self.ax[1].set_ylim(-10, float(self.chisq_ymax_entry.get()))

		self.fig.tight_layout()
		self.canvas.draw()

	def p_of_z_plot(self):

		self.ax[1].clear()
		template_color = '#56B4E9'#'#117733'#'green'
		NIRC_photometry_color = '#D55E00'#'#882255'#'red'
		HST_photometry_color =  '#F786AA'#'#CC6677'#'lightcoral'
		chisq_surface_color = '#E69F00'#'#88CCEE'#'blue'
		zspec_color = '#332288'#'orange'
		alternate_color = 'grey'
		
		object_ID = int(self.id_entry.get())

		bf_output_wavelength, bf_output_flux, bf_output_phot_wavelength, bf_output_phot_template, bf_output_phot_err_template, bf_output_phot, bf_output_phot_err, bf_z_a_value, bf_tempfilt_zgrid, bf_chi2fit = return_eazy_output(object_ID, eazy_self, zout, primary_z = -9999, alt_z = 0.0)
		
		if (self.altz_entry.get() == ''):
			output_wavelength, output_flux, output_phot_wavelength, output_phot_template, output_phot_err_template, output_phot, output_phot_err, z_a_value, tempfilt_zgrid, chi2fit = return_eazy_output(object_ID, eazy_self, zout, primary_z = -9999, alt_z = 0.0)
		else:
			if ((float(self.altz_entry.get()) > eazy_self.param['Z_MIN']) and (float(self.altz_entry.get()) < eazy_self.param['Z_MAX'])):
				output_wavelength, output_flux, output_phot_wavelength, output_phot_template, output_phot_err_template, output_phot, output_phot_err, z_a_value, tempfilt_zgrid, chi2fit = return_eazy_output(object_ID, eazy_self, zout, primary_z = float(self.altz_entry.get()), alt_z = 0.0)
			else:
				self.altz_entry.delete(0, tk.END)
				self.altz_entry.insert(0,"") 
				print("NOTE: Alternate redshift must be between "+str(eazy_self.param['Z_MIN'])+" and "+str(eazy_self.param['Z_MAX']))
				output_wavelength, output_flux, output_phot_wavelength, output_phot_template, output_phot_err_template, output_phot, output_phot_err, z_a_value, tempfilt_zgrid, chi2fit = return_eazy_output(object_ID, eazy_self, zout, primary_z = -9999, alt_z = 0.0)

		p_of_z = np.exp((-1.0) * (chi2fit - np.min(chi2fit))/2.0)
		
		# Now, let's plot the chi-square surface. 
		self.ax[1].plot(tempfilt_zgrid, p_of_z, color = chisq_surface_color, zorder = 40)
		
		# Let's grab the chisq value at the specified redshift. 
		bf_chisq_value = get_chisq(bf_tempfilt_zgrid, bf_chi2fit, bf_z_a_value)
		chisq_value = get_chisq(tempfilt_zgrid, chi2fit, z_a_value)

		objid_index = np.where(eazy_self.OBJID == int(self.current_ID))[0][0]
		z160 = zout['z160'][objid_index]
		z840 = zout['z840'][objid_index]
		
		z500 = zout['z500'][objid_index]

		upper_error = z840 - z500
		lower_error = z500 - z160

		z025 = zout['z025'][objid_index]
		z975 = zout['z975'][objid_index]
		
		self.ax[1].set_xlabel(r'z$_{phot}$')
		self.ax[1].set_ylabel(r'$\chi^2$')

		if (self.altz_entry.get() == ''):
			self.ax[1].axvspan(bf_z_a_value-0.05, bf_z_a_value+0.05, color = template_color, label = 'z = '+str(bf_z_a_value)+', $\chi^2 = $'+str(bf_chisq_value), zorder = 20)
		else:
			self.ax[1].axvspan(bf_z_a_value-0.05, bf_z_a_value+0.05, color = 'grey', label = 'z = '+str(bf_z_a_value)+', $\chi^2 = $'+str(bf_chisq_value), zorder = 19)
			self.ax[1].axvspan(z_a_value-0.05, z_a_value+0.05, color = template_color, label = 'z$_{\mathrm{alt}}$ = '+str(z_a_value)+', $\chi^2 = $'+str(chisq_value), zorder = 20)
		
		self.ax[1].axvspan(z500-0.01, z500+0.01, color = 'red', label = 'z500 = '+str(round(z500,2))+'$^{+'+str(round(upper_error,2))+'}$'+'$_{-'+str(round(lower_error,2))+'}$', alpha = 0.2, zorder = 3)

		self.ax[1].set_xlim(z160-2,z840+2)
		self.ax[1].axvspan(z025, z975, color = 'grey', zorder = 0, alpha = 0.2)
		self.ax[1].axvspan(z160, z840, color = 'grey', zorder = 1, alpha = 0.5)

		pz_legend = self.ax[1].legend(loc = 4, fontsize = 12, frameon=False)
		pz_legend.set_zorder(50)

		if (self.chisq_ymax_entry.get() != ''):
			self.ax[1].set_ylim(-10, float(self.chisq_ymax_entry.get()))

		self.fig.tight_layout()
		self.canvas.draw()
		

  	      
    # This will reset the SED plot axes. 
	def reset_plot(self):

		object_ID = int(self.id_entry.get())

		bf_output_wavelength, bf_output_flux, bf_output_phot_wavelength, bf_output_phot_template, bf_output_phot_err_template, bf_output_phot, bf_output_phot_err, bf_z_a_value, bf_tempfilt_zgrid, bf_chi2fit = return_eazy_output(object_ID, eazy_self, zout, primary_z = -9999, alt_z = 0.0)
		
		if (self.altz_entry.get() == ''):
			output_wavelength, output_flux, output_phot_wavelength, output_phot_template, output_phot_err_template, output_phot, output_phot_err, z_a_value, tempfilt_zgrid, chi2fit = return_eazy_output(object_ID, eazy_self, zout, primary_z = -9999, alt_z = 0.0)
		else:
			output_wavelength, output_flux, output_phot_wavelength, output_phot_template, output_phot_err_template, output_phot, output_phot_err, z_a_value, tempfilt_zgrid, chi2fit = return_eazy_output(object_ID, eazy_self, zout, primary_z = float(self.altz_entry.get()), alt_z = 0.0)

		self.ax[0].set_xlim(0, 5.5)
		self.xmin_entry.delete(0, tk.END)  # Delete current text
		self.xmin_entry.insert(0,0) 
		self.xmax_entry.delete(0, tk.END)  # Delete current text
		self.xmax_entry.insert(0,5.5) 

		pos_fluxes = np.where(output_phot > 0)[0]
		ymin = np.min(output_phot[pos_fluxes]) - (0.4 * np.min(output_phot[pos_fluxes]))
		if (ymin < 0.01):
			ymin = 0.01

		ymax = np.max(output_phot[pos_fluxes]) + (10.0 * np.max(output_phot[pos_fluxes]))
		if (ymax > (100.*np.median(output_phot[pos_fluxes]))):
			ymax = 100.*np.median(output_phot[pos_fluxes])

		self.ax[0].set_ylim(ymin, ymax)

		self.ymin_entry.delete(0, tk.END)  # Delete current text
		self.ymax_entry.delete(0, tk.END)  # Delete current text

		# Now, let's plot the chi-square surface. 
		self.ax[1].plot(tempfilt_zgrid, chi2fit, color = chisq_surface_color, zorder = 40)

		# Let's grab the chisq value at the specified redshift. 
		bf_chisq_value = get_chisq(bf_tempfilt_zgrid, bf_chi2fit, bf_z_a_value)
		chisq_value = get_chisq(tempfilt_zgrid, chi2fit, z_a_value)

		# Got to get the index within the EAZY object
		objid_index = np.where(eazy_self.OBJID == int(self.current_ID))[0][0]
		
		# Here are the output redshift values from the EAZY fit
		z160 = zout['z160'][objid_index]
		z840 = zout['z840'][objid_index]
		
		z500 = zout['z500'][objid_index]

		upper_error = z840 - z500
		lower_error = z500 - z160

		z025 = zout['z025'][objid_index]
		z975 = zout['z975'][objid_index]
		
		self.ax[1].clear()
		self.ax[1].plot(tempfilt_zgrid, chi2fit, color = chisq_surface_color, zorder = 40)
		self.ax[1].set_xlabel(r'z$_{phot}$')
		self.ax[1].set_ylabel(r'$\chi^2$')

		if (self.altz_entry.get() == ''):
			self.ax[1].axvspan(bf_z_a_value-0.05, bf_z_a_value+0.05, color = template_color, label = 'z = '+str(bf_z_a_value)+', $\chi^2 = $'+str(bf_chisq_value), zorder = 20)
		else:
			self.ax[1].axvspan(bf_z_a_value-0.05, bf_z_a_value+0.05, color = 'grey', label = 'z = '+str(bf_z_a_value)+', $\chi^2 = $'+str(bf_chisq_value), zorder = 19)
			self.ax[1].axvspan(z_a_value-0.05, z_a_value+0.05, color = template_color, label = 'z$_{\mathrm{alt}}$ = '+str(z_a_value)+', $\chi^2 = $'+str(chisq_value), zorder = 20)

		self.ax[1].axvspan(z500-0.05, z500+0.05, color = 'red', label = 'z500 = '+str(round(z500,2))+'$^{+'+str(round(upper_error,2))+'}$'+'$_{-'+str(round(lower_error,2))+'}$', alpha = 0.2, zorder = 3)

		self.ax[1].axvspan(z025, z975, color = 'grey', zorder = 0, alpha = 0.1)
		self.ax[1].axvspan(z160, z840, color = 'grey', zorder = 1, alpha = 0.05)

		version_string = return_version_string(args_input_file)

		if (default_convolved == True):
			self.ax[1].set_title('Photometry: '+version_string+', '+str(args_aperture)+', Convolved', fontsize = 15)
		else: 
			self.ax[1].set_title('Photometry: '+version_string+', '+str(args_aperture)+', Not Convolved', fontsize = 15)

		
		chisq_legend = self.ax[1].legend(loc = 2, fontsize = 12, frameon=False)
		chisq_legend.set_zorder(50)

		if (self.chisq_ymax_entry.get() != ''):
			self.ax[1].set_ylim(-10, float(self.chisq_ymax_entry.get()))


		self.fig.tight_layout()
		self.canvas.draw()
    
	def load_images(self):
		self.images = [np.random.rand(50, 50) for _ in range(10)]
		for ax, img in zip(self.image_axes.flatten(), self.images):
			ax.clear()
			ax.imshow(img, cmap='gray')
			ax.axis('off')
		self.fig.tight_layout()
		self.image_canvas.draw()

	def save_eazy_sed(self):
		object_ID = int(self.id_entry.get())

		# Get the photometry version
		version_string = return_version_string(args_input_file)

		if (self.altz_entry.get() == ''):
			output_wavelength, output_flux, output_phot_wavelength, output_phot_template, output_phot_err_template, output_phot, output_phot_err, z_a_value, tempfilt_zgrid, chi2fit = return_eazy_output(object_ID, eazy_self, zout, primary_z = -9999, alt_z = 0.0)
		else:
			output_wavelength, output_flux, output_phot_wavelength, output_phot_template, output_phot_err_template, output_phot, output_phot_err, z_a_value, tempfilt_zgrid, chi2fit = return_eazy_output(object_ID, eazy_self, zout, primary_z = float(self.altz_entry.get()), alt_z = 0.0)

		# Got to update the base final image name with convolved or unconvolved
		output_filename = self.id_entry.get()+"_"+version_string+"_"+str(args_aperture)
		if (default_convolved == True):
			output_filename = output_filename+"_conv"
		else: 
			output_filename = output_filename+"_unconv"

		output_SED_file = output_filename + "_z_"+str(z_a_value) + '_SED_microns_nJy.txt'
		np.savetxt(output_SED_file, np.c_[output_wavelength/1e4, output_flux], fmt = '%f %f', header="Wavelength (microns) Flux (nJy)")	
		print("Saving EAZY SED as "+output_SED_file)

		output_chisq_file = output_filename + '_z_chisq.txt'
		np.savetxt(output_chisq_file, np.c_[tempfilt_zgrid, chi2fit], fmt = '%f %f', header="Redshift chisq")	
		print("Saving EAZY chisq as "+output_chisq_file)


	def save_gui(self):
		# Get the photometry version
		version_string = return_version_string(args_input_file)

		# Need to get the z_a_value to set the base final image name
		if (self.altz_entry.get() == ''):
			output_wavelength, output_flux, output_phot_wavelength, output_phot_template, output_phot_err_template, output_phot, output_phot_err, z_a_value, tempfilt_zgrid, chi2fit = return_eazy_output(int(self.id_entry.get()), eazy_self, zout, primary_z = -9999, alt_z = 0.0)
		else:
			output_wavelength, output_flux, output_phot_wavelength, output_phot_template, output_phot_err_template, output_phot, output_phot_err, z_a_value, tempfilt_zgrid, chi2fit = return_eazy_output(int(self.id_entry.get()), eazy_self, zout, primary_z = float(self.altz_entry.get()), alt_z = 0.0)
		
		# Got to update the base final image name with convolved or unconvolved
		output_filename = self.id_entry.get()+"_"+version_string+"_"+str(args_aperture)
		if (default_convolved == True):
			output_filename = output_filename+"_conv"
		else: 
			output_filename = output_filename+"_unconv"

		# This version has the thumbnails.
		if (args_plot_thumbnails):
			combined_fig, combined_axes = plt.subplots(2, 1, figsize=(13, 9))
			
			# Copy main plot
			self.ax[0].figure.canvas.draw()
			plot_img = np.array(self.ax[0].figure.canvas.renderer.buffer_rgba())
			combined_axes[0].imshow(plot_img)
			combined_axes[0].axis('off')
			
			# Copy image grid
			self.image_fig.canvas.draw()
			image_img = np.array(self.image_fig.canvas.renderer.buffer_rgba())
			combined_axes[1].imshow(image_img)
			combined_axes[1].axis('off')

			# And this final image name also has thumbnails
			output_filename = output_filename +"_z_"+str(z_a_value)+'_SED_thumbnails.png'
			
			combined_fig.tight_layout()
			combined_fig.savefig(output_filename, dpi = 300)
			plt.close(combined_fig)
			print("JADESView Screen saved as "+output_filename)
		else:
			combined_fig, combined_axes = plt.subplots(1, 1, figsize=(11, 4))
			
			# Copy main plot
			self.ax[0].figure.canvas.draw()
			plot_img = np.array(self.ax[0].figure.canvas.renderer.buffer_rgba())
			combined_axes.imshow(plot_img)
			combined_axes.axis('off')
						
			# And this final image name does not indicate that there are thumbnails
			output_filename = output_filename +"_z_"+str(z_a_value)+'_SED.png'
						
			combined_fig.tight_layout()
			combined_fig.savefig(output_filename, dpi = 300)
			plt.close(combined_fig)
			print("JADESView Screen saved as "+output_filename)
		
	def previous_object(self):
		object_iterator = np.where(ID_numbers == int(self.id_entry.get()))[0][0]
		try:
			if (object_iterator > 0):
				object_iterator = object_iterator - 1
				self.id_entry.delete(0, tk.END)  # Delete current text
				self.id_entry.insert(0,str(ID_numbers[object_iterator])) 
				self.update_plot()
			else:
				print("At start of list!")
		except IndexError:
			self.id_entry.delete(0, tk.END)  # Delete current text
			self.id_entry.insert(0,str(ID_numbers[0])) 
			self.update_plot()
			print("At start of list!")
				
	def next_object(self):
		#print(int(self.id_entry.get()), np.where(ID_numbers == int(self.id_entry.get())))
		object_iterator = np.where(ID_numbers == int(self.id_entry.get()))[0][0]
		try:
			#print(object_iterator)
			if (object_iterator < number_objects):
				object_iterator = object_iterator + 1
				self.id_entry.delete(0, tk.END)  # Delete current text
				self.id_entry.insert(0,str(ID_numbers[object_iterator])) 
				self.update_plot()
			else:
				print("At end of list!")
		except IndexError:
			self.id_entry.delete(0, tk.END)  # Delete current text
			self.id_entry.insert(0,str(ID_numbers[-1])) 
			self.update_plot()
			print("At end of list!")
		
	# Saving the individual thumbnails as fits files for each filter shown. 
	def save_thumbnails(self):
		print("Saving thumbnails")

		# Getting the RA/DEC
		object_ID = int(self.id_entry.get())
		object_ID_index = np.where(ID_values == object_ID)[0][0]
		objRA = RA_values[object_ID_index]
		objDEC = DEC_values[object_ID_index]
		output_name_stub = str(object_ID)+'_size_'+str(round(float(self.thumbnail_size.get()),1))+'_arcsec'

		# Creating the output folder, if one doesn't exist 
		output_folder = str(object_ID)+"_Thumbnails"
		if (os.path.exists(output_folder) == False):
			os.mkdir(output_folder)
				
		for q in range(0, number_images):
			image_cutout = return_image_cutout(objRA, objDEC, image_hdu_all[q], image_wcs_all[q], float(self.thumbnail_size.get()))

			hdu_cutout = copy.copy(image_hdu_all[q])
			
			hdu_cutout.data = image_cutout.data
			hdu_cutout.header.update(image_cutout.wcs.to_header())
			hdu_cutout.writeto(output_folder+'/'+output_name_stub+'_'+str(all_images_filter_name[q])+'.fits', overwrite = True)

	def open_fitsmap(self):
		# Getting the RA/DEC
		object_ID = int(self.id_entry.get())
		object_ID_index = np.where(ID_values == object_ID)[0][0]
		objRA = RA_values[object_ID_index]
		objDEC = DEC_values[object_ID_index]

		base_url = fitsmap_link  # Example base URL
		
		# Construct the full URL with query parameters
		full_url = base_url.replace("RA_VALUE", str(objRA))
		full_url = full_url.replace("DEC_VALUE", str(objDEC))
		
		# Open in the default web browser
		webbrowser.open(full_url)

	def on_plot_click(self, event):
		# Check if the click is in the right plot
		if event.inaxes == self.ax[1]: 
			# Ensure click is within the axes
			if event.xdata is not None:  
				# Get X coordinate
				clicked_x = event.xdata 
				# Pass to click handler function
				self.handle_x_click(clicked_x)  

	def handle_x_click(self, x_value):
		# Example action: update a label or print
		self.altz_entry.delete(0, tk.END)
		self.altz_entry.insert(0, f"{x_value:.2f}")
		self.update_plot()
       
	def quit_program(self):
		sys.exit()


######################
#     Arguments      #
######################

parser = argparse.ArgumentParser()

# Redshift list
parser.add_argument(
  '-jv_params','--jv_params',
  help="JADESView input parameter file",
  action="store",
  type=str,
  dest="JADESView_input_file",
  required=True
)
# Redshift list
parser.add_argument(
  '-zlist','--z_number_list',
  help="List of ID numbers and redshfits for a subsample?",
  action="store",
  type=str,
  dest="z_number_list",
  required=False
)

# ID number
parser.add_argument(
  '-id','--id_number',
  help="ID Number?",
  action="store",
  type=int,
  dest="id_number",
  required=False
)

# ID list
parser.add_argument(
  '-idlist','--id_number_list',
  help="List of ID Numbers?",
  action="store",
  type=str,
  dest="id_number_list",
  required=False
)

# command line argument list of objects
parser.add_argument(
  '-idarglist',
  help="Command line argument list of objects",
  action="store",
  type=str,
  dest="idarglist",
  required=False
)


args=parser.parse_args()

if __name__ == '__main__':

	# Read in the various input values from the input file. 
	input_lines = np.loadtxt(args.JADESView_input_file, dtype='str')
	number_input_lines = len(input_lines[:,0])
	for i in range(0, number_input_lines):

		# ANCILLARY FILE FOLDER
		# Ancillary file folder
		if (input_lines[i,0] == 'ancillary_files'):
			ancillary_file_folder = input_lines[i,1]
		
		# CATALOG PARAMETERS
		# input photometry catalog
		if (input_lines[i,0] == 'input_photometry'):
			args_input_file = input_lines[i,1]
		# Convolved photometry? 
		if (input_lines[i,0] == 'convolved'):
			if ((input_lines[i,1] == 'True') or (input_lines[i,1] == 'T') or (input_lines[i,1] == 'Yes') or (input_lines[i,1] == 'Y')):
				default_convolved = True
			else:
				default_convolved = False
		# Which aperture to use?
		if (input_lines[i,0] == 'aperture'):
			args_aperture = input_lines[i,1]
		# Which uncertainty values to use?
		if (input_lines[i,0] == 'uncertainty'):
			default_uncertainty = input_lines[i,1]
		# Which filters to use? 
		if (input_lines[i,0] == 'filters'):
			args_filters = input_lines[i,1]
		# Value of minimum relative error to use in fit?
		if (input_lines[i,0] == 'min_rel_err'):
			args_min_rel_err = float(input_lines[i,1])

		# EAZY PARAMETERS
		# zphot.zeropoint file location
		if (input_lines[i,0] == 'zeropoint'):
			args_zeropoint_filename = input_lines[i,1]
		# template set to use
		if (input_lines[i,0] == 'templates'):
			args_template_param = input_lines[i,1]
		# tempfilt.pkl file (for speeding up the code)
		if (input_lines[i,0] == 'tempfilt'):
			args_tempfilt_filename = input_lines[i,1]
		# Number of cores to use for the EAZY fits
		if (input_lines[i,0] == 'eazy_ncores'):
			args_ezp = int(input_lines[i,1])
		# Flag to use the Asada+24 CGM prescription
		if (input_lines[i,0] == 'asada_cgm'):
			if ((input_lines[i,1] == 'True') or (input_lines[i,1] == 'T') or (input_lines[i,1] == 'Yes') or (input_lines[i,1] == 'Y')):
				args_asada_cgm = True
			else:
				args_asada_cgm = False
		# Fix the redshift to z_spec
		if (input_lines[i,0] == 'fix_zspec'):
			if ((input_lines[i,1] == 'True') or (input_lines[i,1] == 'T') or (input_lines[i,1] == 'Yes') or (input_lines[i,1] == 'Y')):
				fix_zspec = True
			else:
				fix_zspec = False
		
		# input list of image files for creating thumbnails
		if (input_lines[i,0] == 'image_list'):
			all_images_file_name = input_lines[i,1]
		# and, if you want to make thumbnails, setting this flag
		if (input_lines[i,0] == 'plot_thumbnails'):
			if ((input_lines[i,1] == 'True') or (input_lines[i,1] == 'T') or (input_lines[i,1] == 'Yes') or (input_lines[i,1] == 'Y')):
				args_plot_thumbnails = True
			else:
				args_plot_thumbnails = False
		# the size of the thumbnails, in arcseconds
		if (input_lines[i,0] == 'ra_dec_size_value'):
			ra_dec_size_value = float(input_lines[i,1])

		# the raw fitsmap link
		if (input_lines[i,0] == 'fitsmap_link'):
			fitsmap_link = input_lines[i,1]
			fitsmap_link_exists = True

	# CHECKING IF THE FILES EXIST
	# Checking whether the ancillary file folder exists
	if (os.path.isdir(ancillary_file_folder) == False):
		sys.exit("Ancillary file folder doesn't exist: "+ancillary_file_folder)

	# Checking whether the input photometry file exists
	if (os.path.isfile(args_input_file) == False):
		if (os.path.isfile(ancillary_file_folder+'/'+args_input_file) == False):
			sys.exit("Input photometry file doesn't exist: "+args_input_file)

	# Checking whether the filter file exists. It can be linked directly,
	# or placed in the ancillary file folder
	if (os.path.isfile(args_filters) == False):
		if (os.path.isfile(ancillary_file_folder+'/'+args_filters) == False):
			sys.exit("Filters file doesn't exist: "+args_filters)
		else:
			args_filters = ancillary_file_folder+'/'+args_filters

	# Checking whether the template file exists. It can be linked directly,
	# or placed in the ancillary file folder
	if (os.path.isfile(args_template_param) == False):
		if (os.path.isfile(ancillary_file_folder+'/templates/'+args_template_param) == False):
			sys.exit("Template file doesn't exist: "+args_template_param)
		else:
			args_template_param = ancillary_file_folder+'/templates/'+args_template_param

	# The tempfilt filen is not necessary for running the program, but this allows
	# it to be either in the ancillary file folder, or you can point somewhere else.
	if (os.path.isfile(args_tempfilt_filename) == False):
		if (os.path.isfile(ancillary_file_folder+'/'+args_tempfilt_filename) == True):
			args_tempfilt_filename = ancillary_file_folder + '/' + args_tempfilt_filename

	# The zeropoint filename is not necessary for running the program, but this allows
	# it to be either in the ancillary file folder, or you can point somewhere else.
	if (os.path.isfile(args_zeropoint_filename) == False):
		if (os.path.isfile(ancillary_file_folder+'/'+args_zeropoint_filename) == True):
			args_zeropoint_filename = ancillary_file_folder + '/' + args_zeropoint_filename
	
	if (args_plot_thumbnails == True):
		# Checking whether the image list file exists. It can be linked directly,
		# or placed in the ancillary file folder
		if (os.path.isfile(all_images_file_name) == False):
			if (os.path.isfile(ancillary_file_folder+'/'+all_images_file_name) == False):
				sys.exit("Image list file doesn't exist: "+all_images_file_name)
			else:
				all_images_file_name = ancillary_file_folder+'/'+all_images_file_name
	
	
	# If you have a tempfilt.pkl file, you need to import pickle
	if (args_tempfilt_filename):
		import pickle
	
	# The EAZY input file for running the fits		
	EAZY_output_filename = 'EAZY_input_JADESView.dat'


	all_filters_file_name = ancillary_file_folder+'/JADES_All_Filters_EAZYpy.dat'
	filters_all = np.loadtxt(all_filters_file_name, dtype='str')
	all_filters = filters_all[:,0]
	eazy_filter_numbers = filters_all[:,1]
	all_filter_numbers = filters_all[:,1]
	all_header_filters = filters_all[:,2]
	all_filter_waves = filters_all[:,3]
	all_filter_waves = all_filter_waves.astype(float)
	all_filter_bw = filters_all[:,4]
	all_filter_bw = all_filter_bw.astype(float)
	total_all_jades_filters = len(all_filters)

	# Open up the filter file specified by the user. 
	filter_file_name = args_filters
	filters_file = np.loadtxt(filter_file_name, dtype='str')
	filters = filters_file
	number_filters = len(filters)
	
	filter_numbers = np.zeros(number_filters, dtype = int)
	old_eazy_filter_numbers = np.zeros(number_filters, dtype = int)
	header_filters = np.zeros(number_filters, dtype = '|U30')
	short_filters = np.zeros(number_filters, dtype = '|U30')
	filter_waves = np.zeros(number_filters)
	filter_bw = np.zeros(number_filters)
	color_filters = np.zeros(number_filters, dtype = 'U10')
	for n in range(0, number_filters):
		filter_index = np.where(all_filters == filters[n])[0][0]
		header_filters[n] = all_header_filters[filter_index]
		short_filters[n] = all_header_filters[filter_index].replace('HST_','').replace('NRC_','').replace('MIRI_','')
		filter_numbers[n] = all_filter_numbers[filter_index]
		old_eazy_filter_numbers[n] = eazy_filter_numbers[filter_index]
		filter_waves[n] = all_filter_waves[filter_index]
		filter_bw[n] = all_filter_bw[filter_index]
	
		if (header_filters[n].startswith('HST')):
			color_filters[n] = HST_photometry_color
		else:
			color_filters[n] = NIRC_photometry_color
			
	if (args_plot_thumbnails):
		
		# Open up the image list file
		#all_images_file_name = args.image_list
		images_all_txt = np.loadtxt(all_images_file_name, dtype='str')
		if (len(images_all_txt[0]) > 2):
			# First, if the user specifies where the science data extension is, they'll put 
			# them in the second column.
			all_images_filter_name = images_all_txt[:,0]
			all_image_extension_number = images_all_txt[:,1].astype('int')
			all_image_paths = images_all_txt[:,2]
		else:
			# Here, we assume that the science extension is 1 
			all_images_filter_name = images_all_txt[:,0]
			all_image_extension_number = np.zeros(len(all_images_filter_name))+1
			all_image_paths = images_all_txt[:,1]
		
		number_image_filters = len(all_images_filter_name)
		number_images = len(all_image_paths)
		
		image_all = np.empty(0)
		image_hdu_all = np.empty(0)
		image_wcs_all = np.empty(0)
	
		image_short_filters = np.zeros(number_image_filters, dtype = '|U30')
		for i in range(0, number_images):
			image_filter_index = np.where(all_filters == all_images_filter_name[i])[0][0]
	
			image_short_filters[i] = all_header_filters[image_filter_index].replace('HST_','').replace('NRC_','').replace('MIRI_','')
	
			print("Opening up image: "+all_image_paths[i])
			if (all_image_paths[i] == 'NoImage'):
				all_image_paths[i] = 'NoImage.fits'
			image_all = np.append(image_all, fits.open(all_image_paths[i]))
			try:
				image_hdu_all = np.append(image_hdu_all, fits.open(all_image_paths[i])[all_image_extension_number[i]])
			except IndexError:
				print('IndexError')
				image_hdu_all = np.append(image_hdu_all, fits.open(all_image_paths[i]))
				print('Running fits.open('+str(all_image_paths[i])+')')
	
			image_wcs_all = np.append(image_wcs_all, WCS(image_hdu_all[i].header))
			
	print("Opening up Noisy Photometry")
		
	# Now, go from the raw catalog file to the fluxes/errors being explored. 
	ID_values, fluxes_all, errors_all, RA_values, DEC_values = RawCatalog_To_Fluxes_Errs_RA_DEC(args_input_file, short_filters, False, default_convolved, args_aperture, default_uncertainty)
	redshifts = np.zeros(len(ID_values))-9999.0
	number_objects = ID_values.size

	flux_value_cat = fluxes_all
	flux_value_err_cat = errors_all 
	
	if (args_min_rel_err):
		print("Setting minimum relative error.")
		for n in range(0, number_filters):
			for j in range(0, number_objects):
				if ((errors_all[j,n] > 0) & (fluxes_all[j,n] > 0)):
					if ((errors_all[j,n] / fluxes_all[j,n]) < args_min_rel_err):
						errors_all[j,n] = args_min_rel_err * fluxes_all[j,n]

	if (args_plot_thumbnails):
		print("Getting SNR values for the image photometry.")
		image_ID_values, image_flux_value_cat, image_flux_value_err_cat, image_RA_values, image_DEC_values = RawCatalog_To_Fluxes_Errs_RA_DEC(args_input_file, image_short_filters, False, default_convolved, args_aperture, default_uncertainty)
	
		SNR_values = np.zeros([number_objects, number_image_filters])
	
		for j in range(0, number_image_filters):
			
			SNR_values[:,j] = image_flux_value_cat[:,j] / image_flux_value_err_cat[:,j]
			SNR_values[np.argwhere(image_flux_value_err_cat[:,j] < 0),j] = -9999
			SNR_values[np.argwhere(np.isnan(image_flux_value_err_cat[:,j])),j] = -9999

	
	# Create a subsample given an ID number list
	print("Creating ID list subsample.")
	if (args.id_number_list):
		ID_input_file = np.loadtxt(args.id_number_list)
		if (len(ID_input_file.shape) > 1):
			subsample_ID_numbers = ID_input_file[:,0].astype(int)
		else:
			subsample_ID_numbers = ID_input_file.astype(int)
		#ID_numbers = ID_input_file
		number_subsample_objects = len(subsample_ID_numbers)
	
		subsample_id_index = np.zeros(number_subsample_objects, dtype = 'int')
		for j in range(0, number_subsample_objects):
			#print(subsample_ID_numbers[j])
			subsample_id_index[j] = np.where(ID_values == subsample_ID_numbers[j])[0][0]
		
		ID_numbers = ID_values[subsample_id_index]
		#ID_values = ID_values[subsample_id_index]
		redshifts = redshifts[subsample_id_index]
		fluxes_all = fluxes_all[subsample_id_index,:]
		errors_all = errors_all[subsample_id_index,:]
		
		number_objects = len(ID_numbers)
	
	if (args.id_number):
		ID_numbers = np.zeros(1, dtype = int)
		ID_numbers[0] = int(args.id_number)
		#number_input_objects = len(ID_numbers)
		if (args.id_number_list):
			"You can't specify an individual ID and a list, ignoring the list."
	
		subsample_id_index = np.where(ID_values == int(args.id_number))[0][0]
		redshifts = np.array([redshifts[subsample_id_index]])
		fluxes_all = fluxes_all[subsample_id_index,:]
		errors_all = errors_all[subsample_id_index,:]
		
		number_objects = 1
		
	if (args.idarglist):
		ID_numbers = np.array(ast.literal_eval(args.idarglist))
		number_input_objects = len(ID_numbers)

		subsample_id_index = np.zeros(number_input_objects, dtype = 'int')
		for j in range(0, number_input_objects):
			subsample_id_index[j] = np.where(ID_values == ID_numbers[j])[0][0]
		
		#ID_values = ID_values[subsample_id_index]
		redshifts = redshifts[subsample_id_index]
		fluxes_all = fluxes_all[subsample_id_index,:]
		errors_all = errors_all[subsample_id_index,:]
		
		number_objects = len(ID_numbers)

	if not (args.id_number or args.idarglist):
		if not (args.id_number_list):
			if not (args.z_number_list):
				sys.exit("Need to specify an ID number or a list of ID numbers!")

	if (args.z_number_list):
		subsample_ID_numbers = np.loadtxt(args.z_number_list)[:,0]
		redshifts = np.loadtxt(args.z_number_list)[:,1]
		number_subsample_objects = len(subsample_ID_numbers)

		subsample_id_index = np.zeros(number_subsample_objects, dtype = 'int')
		for j in range(0, number_subsample_objects):
			subsample_id_index[j] = np.where(ID_values == subsample_ID_numbers[j])[0][0]
		
		ID_numbers = ID_values[subsample_id_index]
		fluxes_all = fluxes_all[subsample_id_index,:]
		errors_all = errors_all[subsample_id_index,:]
		
		number_objects = len(ID_numbers)

	number_input_objects = len(ID_numbers)	
	
	print('Making EAZY file!')
	
	# Create the EAZY input photometry file. 
	# Write the first two lines of the file
	f = open(EAZY_output_filename, 'w')
	f.write('# id z_spec ')
	for n in range(0, number_filters):
		f.write('f_'+filters[n]+' e_'+filters[n]+' ')
	f.write(' \n')
		
	f.write('# id z_spec ')
	for n in range(0, number_filters):
		f.write('F'+str(filter_numbers[n])+' E'+str(filter_numbers[n])+' ')
	f.write(' \n')
	
	# Write out the data 
	for j in range(0, number_objects):
		f.write(str(int(ID_numbers[j]))+' '+str(redshifts[j])+' ')
		for n in range(0, number_filters):
			if (number_input_objects == 1):
				f.write(str(fluxes_all[n])+' '+str(errors_all[n])+' ')
			else:
				f.write(str(fluxes_all[j,n])+' '+str(errors_all[j,n])+' ')
		f.write(' \n')
	f.close() 
	
	# Write the EAZY translate file (but delete it first)
	exists = os.path.isfile('zphot.translate')
	if exists:
		os.system('rm zphot.translate')
		
	EAZY_translate_file = 'zphot.translate'
	f = open(EAZY_translate_file, 'w')
	for n in range(0, number_filters):
		f.write('f_'+filters[n]+' F'+str(filter_numbers[n])+' \n')
		f.write('e_'+filters[n]+' E'+str(filter_numbers[n])+' \n')
	f.write(' \n')
	f.close() 


	# Set up the parameters
	params = {}
	params['CATALOG_FILE'] = EAZY_output_filename
	params['FILTERS_RES'] = ancillary_file_folder+'/FILTER.RES.latest'
	params['FILTER_FORMAT'] = 1 
	params['SMOOTH_FILTERS'] = False        # Smooth filter curves with Gaussian
	params['SMOOTH_SIGMA'] = 100.         # Gaussian sigma (in Angstroms) to smooth filters

	if (args_template_param):
		params['TEMPLATES_FILE'] = args_template_param # Template definition file

	params['TEMPLATE_COMBOS'] = 'a'     # Template combination options: 
                                        #         1 : one template at a time
                                        #         2 : two templates, read allowed combinations from TEMPLATES_FILE
                                        #        -2 : two templates, all permutations
                                        # a <or> 99 : all templates simultaneously
	params['NMF_TOLERANCE'] = 1.e-4              # Tolerance for non-negative combinations (TEMPLATE_COMBOS=a)
	params['WAVELENGTH_FILE'] = ancillary_file_folder+'/templates/lambda.def' # Wavelength grid definition file
	params['TEMP_ERR_FILE'] = ancillary_file_folder+'/templates/TEMPLATE_ERROR.v2.0.zfourge' # Template error definition file
	params['TEMP_ERR_A2'] = 1.00                # Template error amplitude
	params['SYS_ERR'] = 0.00                    # Systematic flux error (% of flux)
	params['APPLY_IGM'] = True                  # Apply Madau 1995 IGM absorption
	params['ADD_CGM'] = False                    # Add Asada24 CGM damping wing absorption

	if (args_asada_cgm):
		params['ADD_CGM'] = True                    # Add Asada24 CGM damping wing absorption
#		params['SIGMOID_PARAM1'] = 3.4835           # Sigmoid func parameter (A) for the N_HI-z relation in Asada24
#		params['SIGMOID_PARAM2'] = 1.2581           # Sigmoid func parameter (a) for the N_HI-z relation in Asada24
#		params['SIGMOID_PARAM3'] = 18.249           # Sigmoid func parameter (C) for the N_HI-z relation in Asada24
		params['SIGMOID_PARAM1'] = 3.5918           # Sigmoid func parameter (A) for the N_HI-z relation in Asada24 (updated, Feb 20)
		params['SIGMOID_PARAM2'] = 1.8414           # Sigmoid func parameter (a) for the N_HI-z relation in Asada24 (updated, Feb 20)
		params['SIGMOID_PARAM3'] = 18.001           # Sigmoid func parameter (C) for the N_HI-z relation in Asada24 (updated, Feb 20)

	params['SCALE_2175_BUMP'] = 0.00            # Scaling of 2175A bump.  Values 0.13 (0.27) absorb ~10 (20) % at peak.
	params['TEMPLATE_SMOOTH'] = 0.0             # Velocity smoothing (km/s) for templates, < 0 for no smoothing
	params['RESAMPLE_WAVE'] = 'None'

	## Input Files
	params['CATALOG_FORMAT'] = 'ascii.commented_header' # Format if not FITS
	params['MAGNITUDES'] = 'n'                  # Catalog photometry in magnitudes rather than f_nu fluxes
	params['NOT_OBS_THRESHOLD'] = -90           # Ignore flux point if <NOT_OBS_THRESH
	params['N_MIN_COLORS'] = 5                  # Require N_MIN_COLORS to fit
	params['ARRAY_NBITS'] = 32                  # Bit depth of internally-created arrays

	## Output Files
	params['OUTPUT_DIRECTORY'] = 'OUTPUT'       # Directory to put output files in
	params['MAIN_OUTPUT_FILE'] = 'photz'        # Main output file, .zout
	params['PRINT_ERRORS'] = 'y'                # Print 68, 95 and 99% confidence intervals
	params['CHI2_SCALE'] = 1.0                  # Scale ML Chi-squared values to improve confidence intervals
	params['VERBOSE_LOG'] = True                # Dump information from the run into [MAIN_OUTPUT_FILE].param
	params['OBS_SED_FILE'] = True               # Write out observed SED/object, .obs_sed
	params['TEMP_SED_FILE'] = True              # Write out best template fit/object, .temp_sed
	params['POFZ_FILE'] = True                  # Write out Pofz/object, .pz
	params['BINARY_OUTPUT'] = True              # Save OBS_SED, TEMP_SED, PZ in binary format to read with e.g IDL

	## Redshift / Mag prior
	params['APPLY_PRIOR'] = False               # Apply apparent magnitude prior
	params['PRIOR_FILE'] = ancillary_file_folder+'/templates/prior_F160W_TAO.dat' # File containing prior grid
	params['PRIOR_FILTER'] = 448                 # Filter from FILTER_RES corresponding to the columns in PRIOR_FILE, originally 431
	params['PRIOR_ABZP'] = 31.4                 # AB zeropoint of fluxes in catalog.  Needed for calculating apparent mags!
	params['PRIOR_FLOOR'] = 1.e-2

	## Redshift Grid
	if (args.z_number_list):
		params['FIX_ZSPEC'] = True                 # Fix redshift to catalog zspec	
	else:
		params['FIX_ZSPEC'] = False                 # Fix redshift to catalog zspec
	params['Z_MIN'] = 0.00                      # Minimum redshift
	params['Z_MAX'] = 22.0                      # Maximum redshift
	params['Z_STEP'] = 0.01                     # Redshift step size
	params['Z_STEP_TYPE'] = 0                   # 0 = ZSTEP, 1 = Z_STEP*(1+z)

	## Zeropoint Offsets
	params['GET_ZP_OFFSETS'] = False                # Look for zphot.zeropoint file and compute zeropoint offsets
	params['ZP_OFFSET_TOL'] = 1.e-4             # Tolerance for iterative fit for zeropoint offsets [not implemented]

	## Rest-frame colors
	params['REST_FILTERS'] = '---'              # Comma-separated list of rest frame filters to compute
	params['RF_PADDING'] = 1000.                # Padding (Ang) for choosing observed filters around specified rest-frame pair.
	params['RF_ERRORS'] = False                 # Compute RF color errors from p(z)
	params['Z_COLUMN'] = 'z_a'                  # Redshift to use for rest-frame color calculation (z_a, z_p, z_m1, z_m2, z_peak)
	params['USE_ZSPEC_FOR_REST'] = True         # Use z_spec when available for rest-frame colors
	params['READ_ZBIN'] = False                 # Get redshifts from OUTPUT_DIRECTORY/MAIN_OUTPUT_FILE.zbin rather than fitting them.

	## Cosmology
	params['H0'] = 70.0               # Hubble constant (km/s/Mpc)
	params['OMEGA_M'] = 0.3                # Omega_matter
	params['OMEGA_L'] = 0.7                # Omega_lambda

	# Here's where we load in the tempfilt filename, if one exists and is specified. 
	if (args_tempfilt_filename):
		exists = os.path.isfile(args_tempfilt_filename)
		if exists:
			with open(args_tempfilt_filename, 'rb') as inp:
				# Open up the tempfilt file, and compare it to the settings/parameters
				tempfilt_file_data = pickle.load(inp)
				if ((tempfilt_file_data['templates'] == args_template_param) and
					(np.array_equal(tempfilt_file_data['filters'], filters)) and
					(tempfilt_file_data['asada'] == args_asada_cgm) and
					(tempfilt_file_data['Z_MIN'] == params['Z_MIN']) and 
					(tempfilt_file_data['Z_MAX'] == params['Z_MAX']) and
					(tempfilt_file_data['Z_STEP'] == params['Z_STEP']) and
					(tempfilt_file_data['Z_STEP_TYPE'] == params['Z_STEP_TYPE'])
					):
					
					tempfilt_array = tempfilt_file_data['tempfilt_array']
					
				else:
					sys.exit("EXITING: The tempfilt file "+args_tempfilt_filename+" cannot be used in this EAZY configuration")

	if (os.path.isfile(args_zeropoint_filename) == True):
		print("Using zphot.zeropoint file!")
		params['GET_ZP_OFFSETS'] = True                # Look for zphot.zeropoint file and compute zeropoint offsets
		zeropoint_file = args_zeropoint_filename	
		translate_file = 'zphot.translate' 
		#param_file_name = 'zphot.param' 
		if (args_tempfilt_filename):
			exists = os.path.isfile(args_tempfilt_filename)
			if exists:
				eazy_self = eazy.photoz.PhotoZ(params = params, translate_file=translate_file, zeropoint_file=zeropoint_file, load_prior=False, load_products=False, tempfilt = tempfilt_array) 
			else:
				eazy_self = eazy.photoz.PhotoZ(params = params, translate_file=translate_file, zeropoint_file=zeropoint_file, load_prior=False, load_products=False) 
		else:
			eazy_self = eazy.photoz.PhotoZ(params = params, translate_file=translate_file, zeropoint_file=zeropoint_file, load_prior=False, load_products=False) 
	
	else:	
		translate_file = 'zphot.translate' 
		param_file_name = 'zphot.param' 
		if (args_tempfilt_filename):
			exists = os.path.isfile(args_tempfilt_filename)
			if exists:
				eazy_self = eazy.photoz.PhotoZ(params = params, translate_file=translate_file, zeropoint_file=None, load_prior=False, load_products=False, tempfilt = tempfilt_array) 
			else:
				eazy_self = eazy.photoz.PhotoZ(params = params, translate_file=translate_file, zeropoint_file=None, load_prior=False, load_products=False) 
		else:
			eazy_self = eazy.photoz.PhotoZ(params = params, translate_file=translate_file, zeropoint_file=None, load_prior=False, load_products=False) 
	
	if (args_tempfilt_filename):
		exists = os.path.isfile(args_tempfilt_filename)
		if exists:
			print("tempfilt already exists, moving on.")
		else:
			# Let's package up important information.
			tempfilt_file_data = {}
			tempfilt_file_data['templates'] = args_template_param
			tempfilt_file_data['filters'] = filters
			tempfilt_file_data['asada'] = args_asada_cgm
			tempfilt_file_data['Z_MIN'] = params['Z_MIN']
			tempfilt_file_data['Z_MAX'] = params['Z_MAX']
			tempfilt_file_data['Z_STEP'] = params['Z_STEP']
			tempfilt_file_data['Z_STEP_TYPE'] = params['Z_STEP_TYPE']
			tempfilt_file_data['tempfilt_array'] = eazy_self.tempfilt
			with open(args_tempfilt_filename, 'wb') as f:
				#pickle.dump(eazy_self.tempfilt, f)
				pickle.dump(tempfilt_file_data, f)
		
	# Now, we actually fit the catalog, after generating the photometric offsets
	eazy_self.fit_catalog(eazy_self.idx, n_proc=args_ezp)
	
	zout, hdu = eazy_self.standard_output(rf_pad_width=0.5, rf_max_err=2, prior=False, beta_prior=False, save_fits = False)

	# This is the index for the source being shown in the GUI
	object_iterator = 0
	# And here's the ID of the first object in the list, the one that's shown. 
	first_object = ID_numbers[object_iterator]

	
	# And finally, we instantiate the GUI. 
	root = tk.Tk()
	app = PlotGUI(root)
	root.mainloop()
	
	
