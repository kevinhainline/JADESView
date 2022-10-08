#! /usr/bin/env python

import os
import ast
import sys
import math
import time 
import scipy
import argparse
import requests
import tarfile
from requests.auth import HTTPBasicAuth
from io import BytesIO
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table

from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

from astropy.visualization import (MinMaxInterval, ZScaleInterval, LogStretch, ImageNormalize, AsinhStretch, SinhStretch, LinearStretch)

#import matplotlib.backends.tkagg as tkagg
#from matplotlib.backends.backend_agg import FigureCanvasAgg

#from Tkinter import *
#try:
#    from Tkinter import *
#except ImportError:
#    from tkinter import *
#import PIL
#from PIL import ImageTk, Image, ImageGrab

JADESView_input_file = 'JADESView_input_file.dat'

# The default stretch on the various images
defaultstretch = 'LinearStretch'

# The default size of the various images
ra_dec_size_value = 2.0

# The default is to not make the crosshair
make_crosshair = False

def print_nearest_objects(ID_values, catalog_ra, catalog_dec, ra_value, dec_value, distance):
	c = SkyCoord(ra=catalog_ra*u.degree, dec=catalog_dec*u.degree)
	object_ra_dec = SkyCoord(ra=ra_value*u.degree, dec=dec_value*u.degree)
	close_objects = np.where(object_ra_dec.separation(c) < distance*u.arcsec)
	
	print("     ------------------------------------------------------------------------------------------")
	if (len(close_objects[0]) < 1):
		print("     There are no objects within "+str(distance)+" arcseconds of the position.")
	if (len(close_objects[0]) == 1):
		print("      This object is within "+str(distance)+" arcseconds of the position RA = "+str(round(ra_value,6))+", DEC = "+str(round(dec_value,6))+":")
		print("            "+str(ID_values[close_objects[0]][0])+", RA = "+str(catalog_ra[close_objects[0]][0])+", DEC = "+str(catalog_dec[close_objects[0]][0])+", distance = "+str(round(object_ra_dec.separation(c)[close_objects[0][0]].arcsec,3))+" arcsec")
	if (len(close_objects[0]) > 1):
		print("     These objects are within "+str(distance)+" arcseconds of the position RA = "+str(round(ra_value,6))+", DEC = "+str(round(dec_value,6))+":")
		for q in range(0, len(close_objects[0])):
			print("       "+str(ID_values[close_objects[0]][q])+", RA = "+str(catalog_ra[close_objects[0]][q])+", DEC = "+str(catalog_dec[close_objects[0]][q])+", distance = "+str(round(object_ra_dec.separation(c)[close_objects[0][q]].arcsec,3))+" arcsec")
	print("     ------------------------------------------------------------------------------------------")

# 110.7431250 -73.4758056
def parse_ra_dec(ra_dec_string):
	
	ra_dec = ra_dec_string.split(' ')
	
	objRA_list = np.array([-9999.0])
	objDEC_list = np.array([-9999.0])

	if (len(ra_dec) == 6):
		ra_dec_skycoord = SkyCoord(ra_dec[0]+' '+ra_dec[1]+' '+ra_dec[2]+' '
				+ra_dec[3]+' '+ra_dec[4]+' '+ra_dec[5], unit=(u.hourangle, u.deg))
		objRA_list = np.array([ra_dec_skycoord.ra.value])
		objDEC_list = np.array([ra_dec_skycoord.dec.value])
	else:
		if ":" in ra_dec[0]:
			ra_dec[0] = ra_dec[0].replace(":"," ")
		
		if ":" in ra_dec[1]:
			ra_dec[1] = ra_dec[1].replace(":"," ")
		
		if (len(ra_dec[0].split()) == 1):
			objRA_list = np.array([float(ra_dec[0])])
			objDEC_list = np.array([float(ra_dec[1])])
		elif ( (len(ra_dec[0].split()) > 1 ) & (len(ra_dec[0].split()) < 4 )):
			ra_dec_skycoord = SkyCoord(ra_dec[0].split()[0]+' '+ra_dec[0].split()[1]+' '+ra_dec[0].split()[2]+' '
				+ra_dec[1].split()[0]+' '+ra_dec[1].split()[1]+' '+ra_dec[1].split()[2], unit=(u.hourangle, u.deg))
			objRA_list = np.array([ra_dec_skycoord.ra.value])
			objDEC_list = np.array([ra_dec_skycoord.dec.value])
		else:
			objRA_list = np.array([-9999.0])
			objDEC_list = np.array([-9999.0])
		
	return objRA_list, objDEC_list

def findthumbnailradec(ra_value, dec_value, ra_dec_size_value, id_cat, ra_cat, dec_cat):
	#print(ra_value, dec_value, ra_dec_size_value)
	c = SkyCoord(ra=ra_cat*u.degree, dec=dec_cat*u.degree)
	center_of_box = SkyCoord(ra=ra_value*u.degree, dec=dec_value*u.degree)
	d2d = c.separation(center_of_box)
	#idx, d2d, d3d = center_of_box.match_to_catalog_sky(c)

	#214.8405342 52.8179497 10.0

	objects_in_box = np.where(d2d.arcsec < ra_dec_size_value)[0]

	if (len(objects_in_box) > 0):
		return id_cat[objects_in_box], ra_cat[objects_in_box], dec_cat[objects_in_box]
	else:
		return np.array([-9999]), np.array([-9999]), np.array([-9999])

parser = argparse.ArgumentParser()

######################
# Optional Arguments #
######################

# JADESView Input File
parser.add_argument(
  '-input',
  help="JADESView Input File?",
  action="store",
  type=str,
  dest="input",
  required=False
)

# Optional Output Folder
parser.add_argument(
  '-output_folder','--output_folder',
  help="Optional Output Folder?",
  action="store",
  type=str,
  dest="output_folder",
  required=False
)

# Input ra and dec value
parser.add_argument(
  '-radec','--radec',
  help="RA/DEC value (string)",
  action="store",
  type=str,
  dest="radec_value",
  required=False
)

# List of input ra and dec values
parser.add_argument(
  '-radec_list','--radec_list',
  help="RA/DEC value list (string)",
  action="store",
  type=str,
  dest="radec_value_list",
  required=False
)

# Use RA DEC list ID
parser.add_argument(
  '-use_ra_dec_list_id','--use_ra_dec_list_id',
  help="Use ID in RA/DEC list?",
  action="store_true",
  dest="use_ra_dec_list_id",
  required=False
)

# Timer Verbose
parser.add_argument(
  '-create_tarball',
  help="Create tarball for each object?",
  action="store",
  type=str,
  dest="create_tarball",
  required=False
)


args=parser.parse_args()

if (args.input):
	JADESView_input_file = args.input

# Right now, the default canvaswidth is 2000. 
canvaswidth = 2000

# Read in the various input values from the input file. 
input_lines = np.loadtxt(JADESView_input_file, dtype='str')
number_input_lines = len(input_lines[:,0])
for i in range(0, number_input_lines):
	if (input_lines[i,0] == 'input_photometry'):
		input_photometry = input_lines[i,1]
	if (input_lines[i,0] == 'image_list'):
		all_images_file_name = input_lines[i,1]
	if (input_lines[i,0] == 'canvaswidth'):
		canvaswidth = float(input_lines[i,1])
	if (input_lines[i,0] == 'defaultstretch'):
		defaultstretch = input_lines[i,1]
	if (input_lines[i,0] == 'ra_dec_size_value'):
		ra_dec_size_value = float(input_lines[i,1])
	if (input_lines[i,0] == 'fenrir_username'):
		fenrir_username = input_lines[i,1]
	if (input_lines[i,0] == 'fenrir_password'):
		fenrir_password = input_lines[i,1]

# # # # # # # # # # # # # # # # # # 
# Let's open up all the input files

# Open up the photometric catalog
fitsinput = fits.open(input_photometry)
ID_values = fitsinput[1].data['ID'].astype('int')
RA_values = fitsinput[1].data['RA']
DEC_values = fitsinput[1].data['DEC']
number_objects = len(ID_values)

# Open up the image list file
images_all_txt = np.loadtxt(all_images_file_name, dtype='str')
#all_images_filter_name = images_all_txt[:,0]
#all_image_paths = images_all_txt[:,1]
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
#for i in range(0, number_images):
#	#print "Opening up image: "+all_image_paths[i]
#	if (all_image_paths[i] == 'NoImage'):
#		all_image_paths[i] = 'NoImage.fits'
#	image_all = np.append(image_all, fits.open(all_image_paths[i]))
#	try:
#		image_hdu_all = np.append(image_hdu_all, fits.open(all_image_paths[i])[1])
#	except IndexError:
#		image_hdu_all = np.append(image_hdu_all, fits.open(all_image_paths[i]))
#	image_wcs_all = np.append(image_wcs_all, WCS(image_hdu_all[i].header))
for i in range(0, number_images):
	print("Opening up image: "+all_image_paths[i])
	if (all_image_paths[i] == 'NoImage'):
		all_image_paths[i] = 'NoImage.fits'
	image_all = np.append(image_all, fits.open(all_image_paths[i]))
	try:
		#image_hdu_all = np.append(image_hdu_all, fits.open(all_image_paths[i])[1])
		image_hdu_all = np.append(image_hdu_all, fits.open(all_image_paths[i])[all_image_extension_number[i]])
	except IndexError:
		print('IndexError')
		image_hdu_all = np.append(image_hdu_all, fits.open(all_image_paths[i]))
		#print('Running fits.open('+str(all_image_paths[i])+')')
	#print('Running WCS(image_hdu_all[i].header)')
	image_wcs_all = np.append(image_wcs_all, WCS(image_hdu_all[i].header))



sf = canvaswidth / 2000.0 # This is the "shrinkfactor" by which all of the canvas
                          # element positions and sizes are shrunk or expanded. I 
                          # RECOGNIZE THAT I SHOULD PUT THINGS ON A GRID, BUT THAT
                          # WILL COME IN A FUTURE UPDATE, OK

# The fontsize depends on this
fontsize = str(int(20*sf))

canvaswidth = 930

# We use the number of images we have to to set the canvasheight
if (number_images <= 6):
	canvasheight = (canvaswidth*(1.0 / 3.0))  # I lock everything to a 2.35:1 aspect ratio
if ((number_images > 6) & (number_images <= 12)):
	canvasheight = (canvaswidth*(1.0 / 1.9))  # I lock everything to a 2:1 aspect ratio
if (number_images > 12):
	canvasheight = (canvaswidth*(1.0 / 1.5))  # I lock everything to a 1.8:1 aspect ratio

textsizevalue = int(20*sf)
thumbnailsize = 1.5*sf

# We have a few rows of buttons at the bottom of the canvas, and their exact
# positions depend on the number of images. 
if (number_images <= 6):
	toprow_y = 230#720
	bottomrow_y = 300#790
if ((number_images > 6) & (number_images <= 12)):
	toprow_y = 420#870
	bottomrow_y = 490#940
if ((number_images > 12) & (number_images <= 18)):
	toprow_y = 590#1020
	bottomrow_y = 660#1060
if ((number_images > 18) & (number_images <= 32)):
	toprow_y = 740#1020
	bottomrow_y = 810#1060
	thumbnailsize = 1.2*sf

if (args.radec_value):
	ra_dec_string = args.radec_value
	objRA_list, objDEC_list = parse_ra_dec(args.radec_value)
	#ra_dec = ra_dec_string.split(' ')
	#objRA_list = np.array([float(ra_dec[0])])
	#objDEC_list = np.array([float(ra_dec[1])])
	
if (args.radec_value_list):
	if (args.radec_value_list.endswith('.fits')):
		fitsinput = fits.open(args.radec_value_list)
	
		objRA_list = fitsinput[1].data['RA']
		objDEC_list = fitsinput[1].data['DEC']
		
		if (args.use_ra_dec_list_id):
			objID_list = fitsinput[1].data['ID']
		
	else:
		radec_value_list_txt = np.loadtxt(args.radec_value_list, dtype = 'str')

		if (args.use_ra_dec_list_id):
			objID_list = radec_value_list_txt[:,0]
			objRA_list_raw = radec_value_list_txt[:,1]
			objDEC_list_raw = radec_value_list_txt[:,2]

		else:
			objRA_list_raw = radec_value_list_txt[:,0]
			objDEC_list_raw = radec_value_list_txt[:,1]
	
		objRA_list = np.zeros(len(objRA_list_raw))
		objDEC_list = np.zeros(len(objRA_list_raw))
		for i in range(0, len(objRA_list_raw)):
			objRA_list[i], objDEC_list[i] = parse_ra_dec(objRA_list_raw[i]+' '+objDEC_list_raw[i])

number_ra_dec_list = len(objRA_list)
current_ra_dec_index = 0

make_crosshair = False

stretch = defaultstretch

if (args.output_folder):
	output_folder = args.output_folder
else:
	output_folder = './'


for obj in range(0, number_ra_dec_list):

	objRA = objRA_list[obj]
	objDEC = objDEC_list[obj]

	if (args.use_ra_dec_list_id):
		obj_output_file_name = 'Object_'+str(objID_list[obj])
	else:
		obj_output_file_name = 'obj_ra_+'+str(round(objRA,6))+'_dec_'+str(round(objDEC,6))

	print(obj_output_file_name)

	if (not os.path.exists(output_folder+obj_output_file_name+'/')):
		os.makedirs(output_folder+obj_output_file_name+'/')
		#os.makedirs(obj_output_file_name+'/fits/')

	cosdec_center = math.cos(objDEC * 3.141593 / 180.0)
	
	print_nearest_objects(ID_values, RA_values, DEC_values, objRA, objDEC, ra_dec_size_value/2.0)
	
	# Set the position of the object
	position = SkyCoord(str(objRA)+'d '+str(objDEC)+'d', frame='fk5')
	size = u.Quantity((ra_dec_size_value, ra_dec_size_value), u.arcsec)
	
	#fig_photo_objects = np.empty(0, dtype = 'object')

	number_rows = np.floor(number_images/6)+1
	figure_y_min = 2
	#fig, axs = plt.subplots(nrows=number_rows, ncols=6, figsize=(12, 2*number_rows))
	fig = plt.figure(figsize=(12, 2*number_rows))
	
	x_spacer = 0.02
	x_width = (1.0 - (7.0*x_spacer)) / 6.0
	y_spacer = 0.02
	y_width = (1.0 - ((number_rows + 1)*x_spacer)) / number_rows

	for i in range(0, number_images):
		
		# First, let's create the image with all of the thumbnails. 
		#print(all_images_filter_name[i])
		image = image_hdu_all[i].data
		image_hdu = image_hdu_all[i]
		image_wcs = image_wcs_all[i]
		
		#hdu_cutout = image_hdu
		
		row_position = (number_rows + 1) - (np.floor((i)/6)+1)
		col_position = ((i)%6)+1

		ax_pos_x1 = (col_position) * x_spacer + (col_position - 1) * x_width #x_spacer + ((col_position-1) * (x_spacer + x_width))
		ax_pos_x2 = (col_position) * (x_spacer + x_width) #x_spacer + ((col_position-1) * (x_width))
		ax_pos_y1 = (row_position) * y_spacer + (row_position - 1) * y_width #x_spacer + ((col_position-1) * (x_spacer + x_width))
		ax_pos_y2 = (row_position) * (y_spacer + y_width) #x_spacer + ((col_position-1) * (x_width))
				
		#print(i, row_position, col_position, ax_pos_x1, ax_pos_x2, ax_pos_y1, ax_pos_y2, )
		
		# Make the cutout
		#start_time = time.time()
		image_cutout = Cutout2D(image, position, size, wcs=image_wcs)
		#end_time = time.time()
		#print("       Running Cutout2D: " +str(end_time - start_time))

		SNR_fontsize_large = int(15.0*sf)
		SNR_fontsize_small = int(12.0*sf)
						
		# Create the wcs axes
		#plt.clf()
		#fig = plt.figure(figsize=(thumbnailsize,thumbnailsize))
		#ax3 = fig.add_axes([0, 0, 1, 1], projection=image_cutout.wcs)
		ax3 = fig.add_axes([ax_pos_x1, ax_pos_y1, x_width, y_width], projection=image_cutout.wcs)
		
		try:
			ax3.text(0.51, 0.96, all_images_filter_name[i].split('_')[1], transform=ax3.transAxes, fontsize=SNR_fontsize_large, fontweight='bold', ha='center', va='top', color = 'black')
			ax3.text(0.5, 0.95, all_images_filter_name[i].split('_')[1], transform=ax3.transAxes, fontsize=SNR_fontsize_large, fontweight='bold', ha='center', va='top', color = 'white')
		except IndexError:
			ax3.text(0.51, 0.96, all_images_filter_name[i], transform=ax3.transAxes, fontsize=SNR_fontsize_large, fontweight='bold', ha='center', va='top', color = 'black')
			ax3.text(0.5, 0.95, all_images_filter_name[i], transform=ax3.transAxes, fontsize=SNR_fontsize_large, fontweight='bold', ha='center', va='top', color = 'white')
			
		if (make_crosshair == True):
			ax3.plot([0.5, 0.5], [0.65, 0.8], linewidth=2.0, transform=ax3.transAxes, color = 'white')
			ax3.plot([0.5, 0.5], [0.2, 0.35], linewidth=2.0, transform=ax3.transAxes, color = 'white')
			ax3.plot([0.2, 0.35], [0.5, 0.5], linewidth=2.0, transform=ax3.transAxes, color = 'white')
			ax3.plot([0.65, 0.8], [0.5, 0.5], linewidth=2.0, transform=ax3.transAxes, color = 'white')
								
		# Set the color map
		plt.set_cmap('gray')
			
		indexerror = 0
		# Normalize the image using the min-max interval and a square root stretch
		thumbnail = image_cutout.data
		#start_time = time.time()
		if (stretch == 'SinhStretch'):
			try:
				norm = ImageNormalize(thumbnail, interval=ZScaleInterval(), stretch=SinhStretch())
			except IndexError:
				indexerror = 1
			except UnboundLocalError:
				indexerror = 1
		if (stretch == 'LogStretch'):
			try:
				norm = ImageNormalize(thumbnail, interval=ZScaleInterval(), stretch=LogStretch(100))
			except IndexError:
				indexerror = 1
			except UnboundLocalError:
				indexerror = 1
		if (stretch == 'LinearStretch'):
			try:
				norm = ImageNormalize(thumbnail, interval=ZScaleInterval(), stretch=LinearStretch())
			except IndexError:
				indexerror = 1
			except UnboundLocalError:
				indexerror = 1

		
		#start_time = time.time()
		if (all_images_filter_name[i] == 'SEGMAP'):
			ax3.imshow(thumbnail, origin = 'lower', aspect='equal')		
		else:		
			if (indexerror == 0):
				ax3.imshow(thumbnail, origin = 'lower', aspect='equal', norm = norm)
			else:
				ax3.imshow(thumbnail, origin = 'lower', aspect='equal')
		#end_time = time.time()
		#print("       Plotting Thumbnail: " +str(end_time - start_time))
										
		# Keep this handle alive, or else figure will disappear
		#fig_photo_objects = np.append(fig_photo_objects, draw_figure(canvas, fig, loc=(fig_x, fig_y)))
		
		#plt.close('all')
		#end_time = time.time()

		ax3.set_axis_off()

		# Now, let's create the individual thumbnail images. 

		fig2 = plt.figure(figsize=(4, 4))
		ax3 = fig2.add_axes([0, 0, 1, 1], projection=image_cutout.wcs)
		
		try:
			ax3.text(0.51, 0.96, all_images_filter_name[i].split('_')[1], transform=ax3.transAxes, fontsize=SNR_fontsize_large, fontweight='bold', ha='center', va='top', color = 'black')
			ax3.text(0.5, 0.95, all_images_filter_name[i].split('_')[1], transform=ax3.transAxes, fontsize=SNR_fontsize_large, fontweight='bold', ha='center', va='top', color = 'white')
		except IndexError:
			ax3.text(0.51, 0.96, all_images_filter_name[i], transform=ax3.transAxes, fontsize=SNR_fontsize_large, fontweight='bold', ha='center', va='top', color = 'black')
			ax3.text(0.5, 0.95, all_images_filter_name[i], transform=ax3.transAxes, fontsize=SNR_fontsize_large, fontweight='bold', ha='center', va='top', color = 'white')
			
		if (make_crosshair == True):
			ax3.plot([0.5, 0.5], [0.65, 0.8], linewidth=2.0, transform=ax3.transAxes, color = 'white')
			ax3.plot([0.5, 0.5], [0.2, 0.35], linewidth=2.0, transform=ax3.transAxes, color = 'white')
			ax3.plot([0.2, 0.35], [0.5, 0.5], linewidth=2.0, transform=ax3.transAxes, color = 'white')
			ax3.plot([0.65, 0.8], [0.5, 0.5], linewidth=2.0, transform=ax3.transAxes, color = 'white')
								
		# Set the color map
		plt.set_cmap('gray')
			
		indexerror = 0
		# Normalize the image using the min-max interval and a square root stretch
		thumbnail = image_cutout.data
		#start_time = time.time()
		if (stretch == 'SinhStretch'):
			try:
				norm = ImageNormalize(thumbnail, interval=ZScaleInterval(), stretch=SinhStretch())
			except IndexError:
				indexerror = 1
			except UnboundLocalError:
				indexerror = 1
		if (stretch == 'LogStretch'):
			try:
				norm = ImageNormalize(thumbnail, interval=ZScaleInterval(), stretch=LogStretch(100))
			except IndexError:
				indexerror = 1
			except UnboundLocalError:
				indexerror = 1
		if (stretch == 'LinearStretch'):
			try:
				norm = ImageNormalize(thumbnail, interval=ZScaleInterval(), stretch=LinearStretch())
			except IndexError:
				indexerror = 1
			except UnboundLocalError:
				indexerror = 1

		
		#start_time = time.time()
		if (all_images_filter_name[i] == 'SEGMAP'):
			ax3.imshow(thumbnail, origin = 'lower', aspect='equal')		
		else:		
			if (indexerror == 0):
				ax3.imshow(thumbnail, origin = 'lower', aspect='equal', norm = norm)
			else:
				ax3.imshow(thumbnail, origin = 'lower', aspect='equal')

		fig2.savefig(output_folder+obj_output_file_name+'/'+obj_output_file_name+'_'+str(all_images_filter_name[i])+'.png', dpi = 300)
		plt.close(fig2)
		
		# And now, let's save the fits file
		
		#hdu_cutout.data = image_cutout.data
		#hdu_cutout.header.update(image_cutout.wcs.to_header())
		#hdu_cutout.writeto(obj_output_file_name+'/fits/'+obj_output_file_name+'_'+str(all_images_filter_name[i])+'.fits', overwrite = True)

	fig.savefig(output_folder+obj_output_file_name+'/'+obj_output_file_name+'_All_Filters.png', dpi = 300)
	plt.close(fig)
	
	# Create the tarfile
	if (args.create_tarball):
		tar = tarfile.open(output_folder+obj_output_file_name+".tar.gz", "w:gz")
		for name in [output_folder+obj_output_file_name]:
		    tar.add(name)
		tar.close()
		
		# And remove the original file 
