#! /usr/bin/env python

import os
import ast
import sys
import math
import time 
import argparse
import requests
from requests.auth import HTTPBasicAuth
from io import BytesIO
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table

from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

from astropy.visualization import (MinMaxInterval, ZScaleInterval, LogStretch, ImageNormalize, AsinhStretch, SinhStretch, LinearStretch)

import matplotlib.backends.tkagg as tkagg
from matplotlib.backends.backend_agg import FigureCanvasAgg

#from Tkinter import *
try:
    from Tkinter import *
except ImportError:
    from tkinter import *
import PIL
from PIL import ImageTk, Image, ImageGrab

JADESView_input_file = 'JADESView_input_file.dat'

# The default stretch on the various images
defaultstretch = 'LinearStretch'

# The default size of the various images
ra_dec_size_value = 2.0

# The default is to not make the crosshair
make_crosshair = False

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

def resizeimage(image):
	global baseplotwidth
	wpercent = (baseplotwidth / float(image.size[0]))
	hsize = int((float(image.size[1]) * float(wpercent)))
	image = image.resize((baseplotwidth, hsize), PIL.Image.ANTIALIAS)
	photo = ImageTk.PhotoImage(image)
	return photo

def shift_north():
	global fig_photo_objects
	global defaultstretch
	global objRA_list
	global objDEC_list
	global ra_dec_size_value

	height = u.Quantity(ra_dec_size_value, u.arcsec)
	objDEC_list_degree = u.Quantity(objDEC_list[current_ra_dec_index], u.deg)
	objDEC_list_to_shift = objDEC_list_degree + height/5.0
	
	objDEC_list[current_ra_dec_index] = objDEC_list_to_shift.value#objDEC_list[current_ra_dec_index] + 0.00005
	radec_label.configure(text="RA = "+str(np.round(objRA_list[current_ra_dec_index],6))+", DEC = "+str(np.round(objDEC_list[current_ra_dec_index],6)))
	fig_photo_objects = create_thumbnails_ra_dec(canvas, fig_photo_objects, objRA_list[current_ra_dec_index], objDEC_list[current_ra_dec_index], defaultstretch)
	
def shift_south():
	global fig_photo_objects
	global defaultstretch
	global objRA_list
	global objDEC_list
	global ra_dec_size_value

	height = u.Quantity(ra_dec_size_value, u.arcsec)
	objDEC_list_degree = u.Quantity(objDEC_list[current_ra_dec_index], u.deg)
	objDEC_list_to_shift = objDEC_list_degree - height/5.0
	
	objDEC_list[current_ra_dec_index] = objDEC_list_to_shift.value#objDEC_list[current_ra_dec_index] + 0.00005

	radec_label.configure(text="RA = "+str(np.round(objRA_list[current_ra_dec_index],6))+", DEC = "+str(np.round(objDEC_list[current_ra_dec_index],6)))
	fig_photo_objects = create_thumbnails_ra_dec(canvas, fig_photo_objects, objRA_list[current_ra_dec_index], objDEC_list[current_ra_dec_index], defaultstretch)

def shift_east():
	global fig_photo_objects
	global defaultstretch
	global objRA_list
	global objDEC_list
	global ra_dec_size_value

	width = u.Quantity(ra_dec_size_value, u.arcsec)
	objRA_list_degree = u.Quantity(objRA_list[current_ra_dec_index], u.deg)
	objRA_list_to_shift = objRA_list_degree + width/5.0

	objRA_list[current_ra_dec_index] = objRA_list_to_shift.value#objRA_list[current_ra_dec_index] + 0.0002
	radec_label.configure(text="RA = "+str(np.round(objRA_list[current_ra_dec_index],6))+", DEC = "+str(np.round(objDEC_list[current_ra_dec_index],6)))
	fig_photo_objects = create_thumbnails_ra_dec(canvas, fig_photo_objects, objRA_list[current_ra_dec_index], objDEC_list[current_ra_dec_index], defaultstretch)

def shift_west():
	global fig_photo_objects
	global defaultstretch
	global objRA_list
	global objDEC_list
	global ra_dec_size_value

	width = u.Quantity(ra_dec_size_value, u.arcsec)
	objRA_list_degree = u.Quantity(objRA_list[current_ra_dec_index], u.deg)
	objRA_list_to_shift = objRA_list_degree - width/5.0

	objRA_list[current_ra_dec_index] = objRA_list_to_shift.value#objRA_list[current_ra_dec_index] + 0.0002
	radec_label.configure(text="RA = "+str(np.round(objRA_list[current_ra_dec_index],6))+", DEC = "+str(np.round(objDEC_list[current_ra_dec_index],6)))
	fig_photo_objects = create_thumbnails_ra_dec(canvas, fig_photo_objects, objRA_list[current_ra_dec_index], objDEC_list[current_ra_dec_index], defaultstretch)


def linearstretch():
	global sf
	global textsizevalue
	global current_ra_dec_index
	global objRA_list
	global objDEC_list
	global canvas   
	global fig_photo_objects
	global defaultstretch
	global btn5
	global btn6
	global btn7

	btn5.config(height = int(2*sf), width = int(10*sf), fg='black', highlightbackground='white', font=('helvetica bold', textsizevalue))
	btn6.config(height = int(2*sf), width = int(10*sf), fg='grey', highlightbackground='white', font=('helvetica', textsizevalue))
	btn7.config(height = int(2*sf), width = int(10*sf), fg='grey', highlightbackground='white', font=('helvetica', textsizevalue))

	defaultstretch = 'LinearStretch'	
	fig_photo_objects = create_thumbnails_ra_dec(canvas, fig_photo_objects, objRA_list[current_ra_dec_index], objDEC_list[current_ra_dec_index], defaultstretch)

def logstretch():
	global sf
	global textsizevalue
	global current_ra_dec_index
	global objRA_list
	global objDEC_list
	global canvas   
	global fig_photo_objects
	global defaultstretch
	global btn5
	global btn6
	global btn7

	btn5.config(height = int(2*sf), width = int(10*sf), fg='grey', highlightbackground='white', font=('helvetica', textsizevalue))
	btn6.config(height = int(2*sf), width = int(10*sf), fg='black', highlightbackground='white', font=('helvetica bold', textsizevalue))
	btn7.config(height = int(2*sf), width = int(10*sf), fg='grey', highlightbackground='white', font=('helvetica', textsizevalue))

	defaultstretch = 'LogStretch'
	fig_photo_objects = create_thumbnails_ra_dec(canvas, fig_photo_objects, objRA_list[current_ra_dec_index], objDEC_list[current_ra_dec_index], defaultstretch)

def sinhstretch():
	global sf
	global textsizevalue
	global current_ra_dec_index
	global objRA_list
	global objDEC_list
	global canvas   
	global fig_photo_objects
	global defaultstretch
	global btn5
	global btn6
	global btn7

	btn5.config(height = int(2*sf), width = int(10*sf), fg='grey', highlightbackground='white', font=('helvetica', textsizevalue))
	btn6.config(height = int(2*sf), width = int(10*sf), fg='grey', highlightbackground='white', font=('helvetica', textsizevalue))
	btn7.config(height = int(2*sf), width = int(10*sf), fg='black', highlightbackground='white', font=('helvetica bold', textsizevalue))

	defaultstretch = 'SinhStretch'
	fig_photo_objects = create_thumbnails_ra_dec(canvas, fig_photo_objects, objRA_list[current_ra_dec_index], objDEC_list[current_ra_dec_index], defaultstretch)

def changeradecsize():
	global current_ra_dec_index
	global objRA_list
	global objDEC_list
	global canvas   
	global fig_photo_objects
	global ra_dec_size_value
	global defaultstretch
	global e3

	ra_dec_size_value = float(e3.get())
	fig_photo_objects = create_thumbnails_ra_dec(canvas, fig_photo_objects, objRA_list[current_ra_dec_index], objDEC_list[current_ra_dec_index], defaultstretch)

def togglecrosshair():
	global current_ra_dec_index
	global objRA_list
	global objDEC_list
	global canvas   
	global fig_photo_objects
	global ra_dec_size_value
	global defaultstretch

	global make_crosshair
	
	if (make_crosshair == False):
		make_crosshair = True
		btn12.config(font=('helvetica bold', textsizevalue))

	else:
		make_crosshair = False
		btn12.config(font=('helvetica', textsizevalue))
		
	fig_photo_objects = create_thumbnails_ra_dec(canvas, fig_photo_objects, objRA_list[current_ra_dec_index], objDEC_list[current_ra_dec_index], defaultstretch)

def gotoobject():

	global e1
	global current_ra_dec_index
	global objRA_list
	global objDEC_list
	global canvas   
	global fig_photo_objects
	global defaultstretch

	full_ra_dec_string = e1.get()
	objRA_list, objDEC_list = parse_ra_dec(full_ra_dec_string)

	#objRA_list[current_ra_dec_index] = float(full_ra_dec_string.split()[0])
	#objDEC_list[current_ra_dec_index] = float(full_ra_dec_string.split()[1])

	radec_label.configure(text="RA = "+str(np.round(objRA_list[current_ra_dec_index],6))+", DEC = "+str(np.round(objDEC_list[current_ra_dec_index],6)))

	fig_photo_objects = create_thumbnails_ra_dec(canvas, fig_photo_objects, objRA_list[current_ra_dec_index], objDEC_list[current_ra_dec_index], defaultstretch)


def draw_figure(canvas, figure, loc=(0, 0)):
    """ Draw a matplotlib figure onto a Tk canvas

    loc: location of top-left corner of figure on canvas in pixels.
    Inspired by matplotlib source: lib/matplotlib/backends/backend_tkagg.py
    """
    figure_canvas_agg = FigureCanvasAgg(figure)
    figure_canvas_agg.draw()
    figure_x, figure_y, figure_w, figure_h = figure.bbox.bounds
    figure_w, figure_h = int(figure_w), int(figure_h)
    photo = PhotoImage(master=canvas, width=figure_w, height=figure_h)

    # Position: convert from top-left anchor to center anchor
    canvas.create_image(loc[0] + figure_w/2, loc[1] + figure_h/2, image=photo)

    # Unfortunately, there's no accessor for the pointer to the native renderer
    tkagg.blit(photo, figure_canvas_agg.get_renderer()._renderer, colormode=2)

    # Return a handle which contains a reference to the photo object
    # which must be kept live or else the picture disappears
    return photo

def create_thumbnails_ra_dec(canvas, fig_photo_objects, ra_value, dec_value, stretch):
	#global ID_values
	global thumbnailsize
	global ra_dec_size_value
	global image_all
	global image_hdu_all
	global image_wcs_all
	global image_flux_value_err_cat
	global all_images_filter_name
	global number_images
	
	# Let's associate the selected object with it's RA and DEC		
	# Create the object thumbnails. 
	#idx_cat = np.where(ID_values == id_value)[0]
	#idx_cat = id_value_index
	objRA = ra_value
	objDEC = dec_value
			
	cosdec_center = math.cos(objDEC * 3.141593 / 180.0)
		
	# Set the position of the object
	position = SkyCoord(str(objRA)+'d '+str(objDEC)+'d', frame='fk5')
	size = u.Quantity((ra_dec_size_value, ra_dec_size_value), u.arcsec)
	
	fig_photo_objects = np.empty(0, dtype = 'object')
	for i in range(0, number_images):
		image = image_hdu_all[i].data
		image_hdu = image_hdu_all[i]
		image_wcs = image_wcs_all[i]
				
		# Make the cutout
		start_time = time.time()
		image_cutout = Cutout2D(image, position, size, wcs=image_wcs)
		#end_time = time.time()
		#print("       Running Cutout2D: " +str(end_time - start_time))

		SNR_fontsize_large = int(15.0*sf)
		SNR_fontsize_small = int(12.0*sf)
						
		# Create the wcs axes
		plt.clf()
		fig = plt.figure(figsize=(thumbnailsize,thumbnailsize))
		ax3 = fig.add_axes([0, 0, 1, 1], projection=image_cutout.wcs)
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
			
		if (number_images <= 18):
			if (i <= 5):
				fig_x, fig_y = (20*sf)+(175*i*sf), 50*sf
			if ((i > 5) & (i <= 11)):
				fig_x, fig_y = (20*sf)+(175*(i-6)*sf), 225*sf
			if ((i > 11) & (i <= 17)):
				fig_x, fig_y = (20*sf)+(175*(i-12)*sf), 400*sf
		if ((number_images > 18) & (number_images <= 24)):
			if (i <= 7):
				fig_x, fig_y = (20*sf)+(130*i*sf), 50*sf
			if ((i > 7) & (i <= 15)):
				fig_x, fig_y = (20*sf)+(130*(i-8)*sf), 225*sf
			if ((i > 15) & (i <= 23)):
				fig_x, fig_y = (20*sf)+(130*(i-16)*sf), 400*sf
		if ((number_images > 24) & (number_images <= 32)):
			if (i <= 7):
				fig_x, fig_y = (20*sf)+(130*i*sf), 50*sf
			if ((i > 7) & (i <= 15)):
				fig_x, fig_y = (20*sf)+(130*(i-8)*sf), 175*sf
			if ((i > 15) & (i <= 23)):
				fig_x, fig_y = (20*sf)+(130*(i-16)*sf), 300*sf
			if ((i > 23) & (i <= 31)):
				fig_x, fig_y = (20*sf)+(130*(i-24)*sf), 425*sf
							
		# Keep this handle alive, or else figure will disappear
		fig_photo_objects = np.append(fig_photo_objects, draw_figure(canvas, fig, loc=(fig_x, fig_y)))
		plt.close('all')
		end_time = time.time()
#		if (timer_verbose):
#			print("       Plotting Thumbnail: " +str(end_time - start_time))

	return fig_photo_objects

# This is kind of a hack to make sure that the image thumbnail size is printed on
# the output files, since none of the labels or buttons are printed when you use
# canvas.postscript()
def save_canvas():
	global thumbnailsize
	global sf
	global canvas
	global e3

	ra_dec_size_value = float(e3.get())

	fig = plt.figure(figsize=(thumbnailsize*2.51,thumbnailsize/4.0))
	ax8 = fig.add_axes([0, 0, 1, 1])
	ax8.text(0.02, 0.5, "Image Size: "+str(ra_dec_size_value)+"\" x "+str(ra_dec_size_value)+"\"", transform=ax8.transAxes, fontsize=12, fontweight='bold', ha='left', va='center', color = 'black')
	fig_x = 20*sf
	if (number_images <= 6):
		fig_y = 760*sf
	if ((number_images > 6) & (number_images <= 12)):
		fig_y = 900*sf
	if (number_images > 12):
		fig_y = 1040*sf	

	fig_size_object = draw_figure(canvas, fig, loc=(fig_x, fig_y))

	#fig = plt.figure(figsize=(thumbnailsize*2.51,thumbnailsize/4.0))
	fig2 = plt.figure(figsize=(7, 2.5))
	ax9 = fig2.add_axes([0, 0, 1, 1])

	fig2_x = 950
	fig2_y = 650
	
	fig2_size_object = draw_figure(canvas, fig2, loc=(fig2_x, fig2_y))
	#fig.close()

	output_filename = 'Object_RA_'+str(np.round(objRA_list[current_ra_dec_index],6))+'_DEC_'+str(np.round(objDEC_list[current_ra_dec_index],6))+'_JADESView'
	canvas.postscript(file=output_filename+'.eps', colormode='color')
	
	# use PIL to convert to PNG 
	img = Image.open(output_filename+'.eps') 
	os.system('rm '+output_filename+'.eps')
	#print output_filename+'.png'
	img.save(output_filename+'.png', 'png') 

def save_destroy():
	quit()

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

# Timer Verbose
parser.add_argument(
  '-tverb',
  help="Print timer values?",
  action="store",
  type=str,
  dest="tverb",
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

# Open up the image list file
images_all_txt = np.loadtxt(all_images_file_name, dtype='str')
all_images_filter_name = images_all_txt[:,0]
all_image_paths = images_all_txt[:,1]
number_image_filters = len(all_images_filter_name)
number_images = len(all_image_paths)

image_all = np.empty(0)
image_hdu_all = np.empty(0)
image_wcs_all = np.empty(0)
for i in range(0, number_images):
	#print "Opening up image: "+all_image_paths[i]
	if (all_image_paths[i] == 'NoImage'):
		all_image_paths[i] = 'NoImage.fits'
	image_all = np.append(image_all, fits.open(all_image_paths[i]))
	try:
		image_hdu_all = np.append(image_hdu_all, fits.open(all_image_paths[i])[1])
	except IndexError:
		image_hdu_all = np.append(image_hdu_all, fits.open(all_image_paths[i]))
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
	fitsinput = fits.open(args.radec_list)

	objRA_list = fitsinput[1].data['RA']
	objDEC_list = fitsinput[1].data['DEC']

number_ra_dec_list = len(objRA_list)
current_ra_dec_index = 0



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Now, everything is set up, so let's start creating the GUI

# Start by creating the GUI as root
root=Tk()
root.wm_title("JADESView RA/DEC Viewer")

# Create the canvas 
canvas=Canvas(root, height=canvasheight, width=canvaswidth, bg="#ffffff")

canvas.pack(side = TOP, expand=True, fill=BOTH)

# Plot the thumbnails
fig_photo_objects = np.empty(0, dtype = 'object')
fig_photo_objects = create_thumbnails_ra_dec(canvas, fig_photo_objects, objRA_list[current_ra_dec_index], objDEC_list[current_ra_dec_index], defaultstretch)

# # # # # # # # # # # #
# Move to New Object Buttons

radec_label = Label(root, text="RA = "+str(np.round(objRA_list[current_ra_dec_index],6))+", DEC = "+str(np.round(objDEC_list[current_ra_dec_index],6)), fg='black', bg='white', font = "Helvetica "+str(int(textsizevalue*1.5)))
radec_label.place(x=10, y = 9)

if ((args.radec_value_list is not None) & (args.radec_value is None)):
	# Create the Previous Object Button
	btn3 = Button(root, text = 'Previous Object', bd = '5', command = previousobject)  
	btn3.config(height = int(2*sf), width = int(20*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue), padx = 3, pady = 3)
	btn3.place(x = 600*sf, y = bottomrow_y*sf)
	
	# Create the Next Object Button
	btn2 = Button(root, text = 'Next Object', bd = '5', command = nextobject)
	btn2.config(height = int(2*sf), width = int(20*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue), padx = 3, pady = 3)
	btn2.place(x = 917*sf, y = bottomrow_y*sf)

if ((args.radec_value_list is None) & (args.radec_value is not None)):
	# Create the Object Entry Field and Button
	ra_dec_entry=StringVar()
	Label(root, text="RA/DEC: ", font=('helvetica', textsizevalue), fg="#000000", bg="#ffffff").place(x=600*sf, y = (bottomrow_y-25.0)*sf)
	e1 = Entry(root, width = int(35*sf), textvariable = ra_dec_entry, font=('helvetica', textsizevalue), fg="#000000", bg="#ffffff")
	e1.place(x = 700*sf, y = (bottomrow_y-30.0)*sf)

	btn9 = Button(root, text = 'Go', bd = '5', command = gotoobject)  
	btn9.config(height = int(1*sf), width = int(4*sf), fg='blue', highlightbackground='white', font=('helvetica', textsizevalue))
	#btn2.pack(side = 'bottom')
	btn9.place(x = 600*sf, y = (bottomrow_y+9.0)*sf)

# # # # # # # # # # # # # # # # # 
# Move Up/Down/Left/Right Button

btn13 = Button(root, text = 'N', bd = '5', command = shift_north)  
btn13.config(height = int(2*sf), width = int(2*sf), fg='red', highlightbackground='white', font=('helvetica', textsizevalue))
btn13.place(x = 600*sf, y = (toprow_y-30.0)*sf)

btn14 = Button(root, text = 'S', bd = '5', command = shift_south)  
btn14.config(height = int(2*sf), width = int(2*sf), fg='red', highlightbackground='white', font=('helvetica', textsizevalue))
btn14.place(x = 600*sf, y = (toprow_y+10.0)*sf)

btn15 = Button(root, text = 'E', bd = '5', command = shift_east)  
btn15.config(height = int(2*sf), width = int(2*sf), fg='red', highlightbackground='white', font=('helvetica', textsizevalue))
btn15.place(x = 540*sf, y = (toprow_y-12.0)*sf)

btn16 = Button(root, text = 'W', bd = '5', command = shift_west)  
btn16.config(height = int(2*sf), width = int(2*sf), fg='red', highlightbackground='white', font=('helvetica', textsizevalue))
btn16.place(x = 660*sf, y = (toprow_y-12.0)*sf)



# # # # # # # # # # # #
# Quit Button

btn4 = Button(root, text = 'Quit', bd = '5', command = save_destroy)  
btn4.config(height = int(2*sf), width = int(10*sf), fg='red', highlightbackground='white', font=('helvetica', textsizevalue))
btn4.place(x = 900*sf, y = (bottomrow_y+10.0)*sf)

# # # # # # # # # # # #
# Save Canvas Button

# save_canvas_imagegrab()
btn4 = Button(root, text = 'Save Canvas', bd = '5', command = save_canvas)  
btn4.config(height = int(2*sf), width = int(15*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue))
btn4.place(x = 700*sf, y = (bottomrow_y+10.0)*sf)

# # # # # # # # # # # #
# Image Stretch Buttons

Label(root, text="Stretch", font = ('helvetica', int(20*sf)), fg="#000000", bg='#ffffff').place(x=20*sf, y = (bottomrow_y+ 14.0)*sf)

# Create the LinearStretch Button
btn5 = Button(root, text = 'Linear', bd = '5', command = linearstretch)  
if (defaultstretch == 'LinearStretch'):
	btn5.config(height = int(2*sf), width = int(10*sf), fg='black', highlightbackground='white', font=('helvetica bold', textsizevalue))
else:
	btn5.config(height = int(2*sf), width = int(10*sf), fg='grey', highlightbackground='white', font=('helvetica', textsizevalue))
btn5.place(x = 100*sf, y = (bottomrow_y+10.0)*sf)


# Create the LogStretch Button
btn6 = Button(root, text = 'Log', bd = '5', command = logstretch)  
if (defaultstretch == 'LogStretch'):
	btn6.config(height = int(2*sf), width = int(10*sf), fg='black', highlightbackground='white', font=('helvetica bold', textsizevalue))
else:
	btn6.config(height = int(2*sf), width = int(10*sf), fg='grey', highlightbackground='white', font=('helvetica', textsizevalue))
btn6.place(x = 250*sf, y = (bottomrow_y+10.0)*sf)


# Create the Sinh Button
btn7 = Button(root, text = 'Sinh', bd = '5', command = sinhstretch)  
if (defaultstretch == 'SinhStretch'):
	btn7.config(height = int(2*sf), width = int(10*sf), fg='black', highlightbackground='white', font=('helvetica bold', textsizevalue))
else:
	btn7.config(height = int(2*sf), width = int(10*sf), fg='grey', highlightbackground='white', font=('helvetica', textsizevalue))
btn7.place(x = 400*sf, y = (bottomrow_y+10.0)*sf)


# Create the RA and DEC size field 
Label(root, text="Box size", font=('helvetica', textsizevalue), fg="#000000", bg="#ffffff").place(x=20*sf, y = (toprow_y+40.0)*sf)
e3 = Entry(root, width = int(7*sf), font=('helvetica', textsizevalue), fg="#000000", bg="#ffffff")
e3.place(x = 120*sf, y = (toprow_y+40.0)*sf)
e3.insert(0, str(ra_dec_size_value))
Label(root, text="arcsec", font=('helvetica', textsizevalue), fg="#000000", bg="#ffffff").place(x=210*sf, y = (toprow_y+42.0)*sf)
btn8 = Button(root, text = 'Resize', bd = '5', command = changeradecsize)  
btn8.config(height = 1, width = int(5*sf), fg='blue', highlightbackground = 'white', font=('helvetica', textsizevalue))
btn8.place(x = 280*sf, y = (toprow_y+37.0)*sf)

btn12 = Button(root, text = 'Crosshair', bd = '5', command = togglecrosshair)  
btn12.config(height = 1, width = int(10*sf), fg='blue', highlightbackground = 'white', font=('helvetica', textsizevalue))
btn12.place(x = 400*sf, y = (toprow_y+37.0)*sf)

root.mainloop()