import os
import ast
import sys
import math
import argparse
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
from PIL import ImageTk, Image

JADESView_input_file = 'JADESView_input_file.dat'

# The default stretch on the various images
defaultstretch = 'LogStretch'

# The default size of the various images
ra_dec_size_value = 2.0


def resizeimage(image):
	global baseplotwidth
	wpercent = (baseplotwidth / float(image.size[0]))
	hsize = int((float(image.size[1]) * float(wpercent)))
	image = image.resize((baseplotwidth, hsize), PIL.Image.ANTIALIAS)
	photo = ImageTk.PhotoImage(image)
	return photo
	

def highz():
	global current_index
	global ID_iterator
	global ID_list
	global highZflag_array
	
	highZflag_array[current_index] = 1
	current_id = ID_list[ID_iterator]
	print("Object "+str(current_id)+" is a high-redshift candidate.")

def badfit():
	global current_index
	global ID_iterator
	global ID_list
	global badfitflag_array
	
	badfitflag_array[current_index] = 1
	current_id = ID_list[ID_iterator]
	print("Object "+str(current_id)+" has a bad fit.")

def baddata():
	global current_index
	global ID_iterator
	global ID_list
	global baddataflag_array
	
	baddataflag_array[current_index] = 1
	current_id = ID_list[ID_iterator]
	print("Object "+str(current_id)+" object has bad data.")


def nextobject():
	global e2
	global ID_iterator
	global current_index
	global ID_list
	global ID_list_indices
	global photo
	global new_photo
	global item4
	global item5
	global canvas   
	global fig_photo_objects
	global defaultstretch

	global eazy_positionx, eazy_positiony
	global eazytext_positionx, eazytext_positiony
	global beagle_positionx, beagle_positiony
	global beagletext_positionx, beagletext_positiony

	notes_values[current_index] = e2.get()
	e2.delete(0,END)

	canvas.delete(item4)
	canvas.delete(item5)
	if (ID_iterator < len(ID_list)-1):
		ID_iterator = ID_iterator+1
#	else:
#		len(ID_list)-1
#		save_destroy()
	
	current_index = ID_list_indices[ID_iterator]
	current_id = ID_list[ID_iterator]
	e2.insert(0, notes_values[current_index])

	image = Image.open(EAZY_files+str(current_id)+"_EAZY_SED.png")
	image = cropEAZY(image)
	photo = resizeimage(image)
	item4 = canvas.create_image(eazy_positionx, eazy_positiony, image=photo)
	
	new_image = Image.open(BEAGLE_files+str(current_id)+"_BEAGLE_SED.png")
	new_photo = resizeimage(new_image)
	item5 = canvas.create_image(beagle_positionx, beagle_positiony, image=new_photo)

	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, current_id, current_index, defaultstretch)

def previousobject():
	global ID_iterator
	global current_index
	global e2
	global ID_list
	global ID_list_indices
	global photo
	global new_photo
	global item4
	global item5
	global canvas   
	global fig_photo_objects
	global defaultstretch

	global eazy_positionx, eazy_positiony
	global eazytext_positionx, eazytext_positiony
	global beagle_positionx, beagle_positiony
	global beagletext_positionx, beagletext_positiony

	notes_values[current_index] = e2.get()
	e2.delete(0,END)

	canvas.delete(item4)
	canvas.delete(item5)
	if (ID_iterator > 0):
		ID_iterator = ID_iterator-1
	else:
		ID_iterator = 0
	
	current_index = ID_list_indices[ID_iterator]
	current_id = ID_list[ID_iterator]
	e2.insert(0, notes_values[current_index])

	image = Image.open(EAZY_files+str(current_id)+"_EAZY_SED.png")
	image = cropEAZY(image)
	photo = resizeimage(image)
	item4 = canvas.create_image(eazy_positionx, eazy_positiony, image=photo)
	
	new_image = Image.open(BEAGLE_files+str(current_id)+"_BEAGLE_SED.png")
	new_photo = resizeimage(new_image)
	item5 = canvas.create_image(beagle_positionx, beagle_positiony, image=new_photo)

	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, current_id, current_index, defaultstretch)


def gotoobject():

	global ID_iterator
	global current_index
	global e2
	global ID_list
	global ID_list_indices
	global photo
	global new_photo
	global item4
	global item5
	global canvas   
	global fig_photo_objects
	global defaultstretch

	global eazy_positionx, eazy_positiony
	global eazytext_positionx, eazytext_positiony
	global beagle_positionx, beagle_positiony
	global beagletext_positionx, beagletext_positiony

	notes_values[current_index] = e2.get()
	e2.delete(0,END)

	if (e1.get().isdigit() == True):
		canvas.delete(item4)
		canvas.delete(item5)

		#current_id = int(e1.get())
		ID_iterator = np.where(ID_list == int(e1.get()))[0][0]

		current_index = ID_list_indices[ID_iterator]
		current_id = ID_list[ID_iterator]
		e2.insert(0, notes_values[current_index])
	
		image = Image.open(EAZY_files+str(current_id)+"_EAZY_SED.png")
		image = cropEAZY(image)
		photo = resizeimage(image)
		item4 = canvas.create_image(eazy_positionx, eazy_positiony, image=photo)
		
		new_image = Image.open(BEAGLE_files+str(current_id)+"_BEAGLE_SED.png")
		new_photo = resizeimage(new_image)
		item5 = canvas.create_image(beagle_positionx, beagle_positiony, image=new_photo)
	
		fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, current_id, current_index, defaultstretch)
	else:
		print("That's not a valid ID number.")


# This will remove the thumbnails, for future work
def cropEAZY(img):

	#output_image = img.crop((0, 0, 3300, 1600))
	output_image = img.crop((0, 0, 3300, 1480))

	return output_image

def linearstretch():
	global sf
	global textsizevalue
	global ID_iterator
	global ID_list
	global ID_list_indices
	global canvas   
	global fig_photo_objects
	global defaultstretch
	global btn5
	global btn6
	global btn7

	btn5.config(height = int(2*sf), width = int(10*sf), fg='black', font=('helvetica', textsizevalue))
	btn6.config(height = int(2*sf), width = int(10*sf), fg='grey', font=('helvetica', textsizevalue))
	btn7.config(height = int(2*sf), width = int(10*sf), fg='grey', font=('helvetica', textsizevalue))

	defaultstretch = 'LinearStretch'	
	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, ID_list[ID_iterator], ID_list_indices[ID_iterator], defaultstretch)

def logstretch():
	global sf
	global textsizevalue
	global ID_iterator
	global ID_list
	global ID_list_indices
	global canvas   
	global fig_photo_objects
	global defaultstretch
	global btn5
	global btn6
	global btn7

	btn5.config(height = int(2*sf), width = int(10*sf), fg='grey', font=('helvetica', textsizevalue))
	btn6.config(height = int(2*sf), width = int(10*sf), fg='black', font=('helvetica', textsizevalue))
	btn7.config(height = int(2*sf), width = int(10*sf), fg='grey', font=('helvetica', textsizevalue))

	defaultstretch = 'LogStretch'
	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, ID_list[ID_iterator], ID_list_indices[ID_iterator], defaultstretch)

def asinhstretch():
	global sf
	global textsizevalue
	global ID_iterator
	global ID_list
	global ID_list_indices
	global canvas   
	global fig_photo_objects
	global defaultstretch
	global btn5
	global btn6
	global btn7

	btn5.config(height = int(2*sf), width = int(10*sf), fg='grey', font=('helvetica', textsizevalue))
	btn6.config(height = int(2*sf), width = int(10*sf), fg='grey', font=('helvetica', textsizevalue))
	btn7.config(height = int(2*sf), width = int(10*sf), fg='black', font=('helvetica', textsizevalue))

	defaultstretch = 'AsinhStretch'
	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, ID_list[ID_iterator], ID_list_indices[ID_iterator], defaultstretch)

# This is kind of a hack to make sure that the image thumbnail size is printed on
# the output files, since none of the labels or buttons are printed when you use
# canvas.postscript()
def save_canvas():
	global thumbnailsize
	global sf
	global canvas
	global ID_list
	global ID_iterator
	global e3

	ra_dec_size_value = float(e3.get())

	fig = plt.figure(figsize=(thumbnailsize*2.51,thumbnailsize/4.0))
	ax8 = fig.add_axes([0, 0, 1, 1])
	ax8.text(0.1, 0.5, "Image Size: "+str(ra_dec_size_value)+"\" x "+str(ra_dec_size_value)+"\"", transform=ax8.transAxes, fontsize=12, fontweight='bold', ha='left', va='center', color = 'black')
	
	fig_x = 20*sf
	if (number_images <= 6):
		fig_y = 760*sf
	if ((number_images > 6) & (number_images <= 12)):
		fig_y = 900*sf
	if (number_images > 12):
		fig_y = 1040*sf	

	fig_size_object = draw_figure(canvas, fig, loc=(fig_x, fig_y))

	current_id = ID_list[ID_iterator]

	output_filename = str(current_id)+'_JADESView'
	canvas.postscript(file=output_filename+'.eps', colormode='color')
	# use PIL to convert to PNG 
	img = Image.open(output_filename+'.eps') 
	#os.system('rm '+output_filename+'.eps')
	#print output_filename+'.png'
	#img.save(output_filename+'.png', 'png') 


def changeradecsize():
	global ID_iterator
	global ID_list
	global ID_list_indices
	global canvas   
	global fig_photo_objects
	global ra_dec_size_value
	global defaultstretch
	global e3

	ra_dec_size_value = float(e3.get())
	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, ID_list[ID_iterator], ID_list_indices[ID_iterator], defaultstretch)


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

def create_thumbnails(canvas, fig_photo_objects, id_value, id_value_index, stretch):
	#global ID_values
	global thumbnailsize
	global RA_values
	global DEC_values
	global ra_dec_size_value
	global image_all
	global image_hdu_all
	global image_wcs_all
	global image_flux_value_err_cat
	global all_images_filter_name
	global number_images
	global SNR_values
	
	# Let's associate the selected object with it's RA and DEC		
	# Create the object thumbnails. 
	#idx_cat = np.where(ID_values == id_value)[0]
	idx_cat = id_value_index
	objRA = RA_values[idx_cat]
	objDEC = DEC_values[idx_cat]
			
	cosdec_center = math.cos(objDEC * 3.141593 / 180.0)
		
	# Set the position of the object
	position = SkyCoord(str(objRA)+'d '+str(objDEC)+'d', frame='fk5')
	size = u.Quantity((ra_dec_size_value, ra_dec_size_value), u.arcsec)
	
	fig_photo_objects = np.empty(0, dtype = 'object')
	for i in range(0, number_images):
		image = image_all[i].data
		image_hdu = image_hdu_all[i]
		image_wcs = image_wcs_all[i]
				
		if (all_images_filter_name[i] == 'HST_F814W'):
			image_wcs.sip = None
		
		if (image_flux_value_err_cat[idx_cat, i] > -9999):
			# Make the cutout
			image_cutout = Cutout2D(image, position, size, wcs=image_wcs)
			
			# Create the wcs axes
			plt.clf()
			fig = plt.figure(figsize=(thumbnailsize,thumbnailsize))
			ax3 = fig.add_axes([0, 0, 1, 1], projection=image_cutout.wcs)
			ax3.text(0.51, 0.96, all_images_filter_name[i].split('_')[1], transform=ax3.transAxes, fontsize=12, fontweight='bold', ha='center', va='top', color = 'black')
			ax3.text(0.5, 0.95, all_images_filter_name[i].split('_')[1], transform=ax3.transAxes, fontsize=12, fontweight='bold', ha='center', va='top', color = 'white')
			if (number_images <= 18):
				if (SNR_values[idx_cat, i] > -100):
					ax3.text(0.96, 0.06, 'SNR = '+str(round(SNR_values[idx_cat, i],2)), transform=ax3.transAxes, fontsize=12, fontweight='bold', horizontalalignment='right', color = 'black')
					ax3.text(0.95, 0.05, 'SNR = '+str(round(SNR_values[idx_cat, i],2)), transform=ax3.transAxes, fontsize=12, fontweight='bold', horizontalalignment='right', color = 'white')
				else:
					ax3.text(0.96, 0.06, 'SNR < -100', transform=ax3.transAxes, fontsize=12, fontweight='bold', horizontalalignment='right', color = 'black')
					ax3.text(0.95, 0.05, 'SNR < -100', transform=ax3.transAxes, fontsize=12, fontweight='bold', horizontalalignment='right', color = 'white')
			else:
				if (SNR_values[idx_cat, i] > -100):
					ax3.text(0.96, 0.06, 'SNR = '+str(round(SNR_values[idx_cat, i],2)), transform=ax3.transAxes, fontsize=9, fontweight='bold', horizontalalignment='right', color = 'black')
					ax3.text(0.95, 0.05, 'SNR = '+str(round(SNR_values[idx_cat, i],2)), transform=ax3.transAxes, fontsize=9, fontweight='bold', horizontalalignment='right', color = 'white')
				else:
					ax3.text(0.96, 0.06, 'SNR < -100', transform=ax3.transAxes, fontsize=9, fontweight='bold', horizontalalignment='right', color = 'black')
					ax3.text(0.95, 0.05, 'SNR < -100', transform=ax3.transAxes, fontsize=9, fontweight='bold', horizontalalignment='right', color = 'white')
						
			# Set the color map
			plt.set_cmap('gray')
			
			indexerror = 0		
			# Normalize the image using the min-max interval and a square root stretch
			thumbnail = image_cutout.data
			if (stretch == 'AsinhStretch'):
				try:
					norm = ImageNormalize(thumbnail, interval=ZScaleInterval(), stretch=AsinhStretch())
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
			
			if (indexerror == 0):
				ax3.imshow(thumbnail, origin = 'lower', aspect='equal', norm = norm)
			else:
				ax3.imshow(thumbnail, origin = 'lower', aspect='equal')
			
			if (number_images <= 18):
				if (i <= 5):
					fig_x, fig_y = (20*sf)+(175*i*sf), 500*sf
				if ((i > 5) & (i <= 11)):
					fig_x, fig_y = (20*sf)+(175*(i-6)*sf), 675*sf
				if ((i > 11) & (i <= 17)):
					fig_x, fig_y = (20*sf)+(175*(i-12)*sf), 850*sf
			if ((number_images > 18) & (number_images <= 24)):
				if (i <= 7):
					fig_x, fig_y = (20*sf)+(130*i*sf), 500*sf
				if ((i > 7) & (i <= 15)):
					fig_x, fig_y = (20*sf)+(130*(i-8)*sf), 675*sf
				if ((i > 15) & (i <= 23)):
					fig_x, fig_y = (20*sf)+(130*(i-16)*sf), 850*sf
			if ((number_images > 24) & (number_images <= 32)):
				if (i <= 7):
					fig_x, fig_y = (20*sf)+(130*i*sf), 500*sf
				if ((i > 7) & (i <= 15)):
					fig_x, fig_y = (20*sf)+(130*(i-8)*sf), 625*sf
				if ((i > 15) & (i <= 23)):
					fig_x, fig_y = (20*sf)+(130*(i-16)*sf), 750*sf
				if ((i > 23) & (i <= 31)):
					fig_x, fig_y = (20*sf)+(130*(i-24)*sf), 875*sf
							
			# Keep this handle alive, or else figure will disappear
			fig_photo_objects = np.append(fig_photo_objects, draw_figure(canvas, fig, loc=(fig_x, fig_y)))
			plt.close('all')

	return fig_photo_objects



def save_destroy():
	global ID_values
	global current_index
	global ID_iterator
	global highZflag_array
	global output_flags_file
	global output_notes_file
	global e2 
	
	if (os.path.exists(output_flags_file)):
		os.system('rm '+output_flags_file)
	if (os.path.exists(output_notes_file)):
		os.system('rm '+output_notes_file)

	# First, let's make the dtype and colnames arrays
	colnames = np.zeros(4, dtype ='S20')
	dtype = np.zeros(4, dtype ='str')
	colnames[0] = 'ID'
	colnames[1] = 'HighZFlag'
	colnames[2] = 'BadFitFlag'
	colnames[3] = 'BadDataFlag'
	dtype[0] = 'I'
	dtype[1] = 'I'
	dtype[2] = 'I'
	dtype[3] = 'I'
	
	# And now let's assemble the data array
	output_data = np.zeros([number_objects, 4])
	output_data[:,0] = ID_values
	output_data[:,1] = highZflag_array
	output_data[:,2] = badfitflag_array
	output_data[:,3] = baddataflag_array
	
	# And finally, let's write out the output file.
	outtab = Table(output_data, names=colnames, dtype=dtype)
	outtab.write(output_flags_file)

	#notes_values[ID_iterator] = e2.get()
	current_id = ID_list[ID_iterator]
	notes_values[current_index] = e2.get()

	w = open(output_notes_file, 'a')
	w.write('#ID    Notes \n')
	for z in range(0, len(ID_values)):
		if (notes_values[z] != ''):
			w.write(str(ID_values[z])+'    '+str(notes_values[z])+'\n')
	w.close()

	quit()
	#root.destroy()



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


# ID number
parser.add_argument(
  '-id',
  help="ID Number?",
  action="store",
  type=int,
  dest="id_number",
  required=False
)

# ID list
parser.add_argument(
  '-idlist',
  help="List of ID Numbers?",
  action="store",
  type=str,
  dest="id_number_list",
  required=False
)

# User Depths 
parser.add_argument(
  '-idarglist',
  help="Command line argument list of objects",
  action="store",
  type=str,
  dest="idarglist",
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
	if (input_lines[i,0] == 'EAZY_files'):
		EAZY_files = input_lines[i,1]
	if (input_lines[i,0] == 'BEAGLE_files'):
		BEAGLE_files = input_lines[i,1]
	if (input_lines[i,0] == 'output_flags_file'):
		output_flags_file = input_lines[i,1]
	if (input_lines[i,0] == 'output_notes_file'):
		output_notes_file = input_lines[i,1]
	if (input_lines[i,0] == 'canvaswidth'):
		canvaswidth = float(input_lines[i,1])
	if (input_lines[i,0] == 'defaultstretch'):
		defaultstretch = input_lines[i,1]
	if (input_lines[i,0] == 'ra_dec_size_value'):
		ra_dec_size_value = float(input_lines[i,1])
	

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
	image_hdu_all = np.append(image_hdu_all, fits.open(all_image_paths[i])[0])
	image_wcs_all = np.append(image_wcs_all, WCS(image_hdu_all[i].header))


sf = canvaswidth / 2000.0 # This is the "shrinkfactor" by which all of the canvas
                          # element positions and sizes are shrunk or expanded. I 
                          # RECOGNIZE THAT I SHOULD PUT THINGS ON A GRID, BUT THAT
                          # WILL COME IN A FUTURE UPDATE, OK

# We use the number of images we have to to set the canvasheight
if (number_images <= 6):
	canvasheight = (canvaswidth*(1.0 / 2.35))  # I lock everything to a 2.35:1 aspect ratio
if ((number_images > 6) & (number_images <= 12)):
	canvasheight = (canvaswidth*(1.0 / 2.0))  # I lock everything to a 2:1 aspect ratio
if (number_images > 12):
	canvasheight = (canvaswidth*(1.0 / 1.8))  # I lock everything to a 1.8:1 aspect ratio

baseplotwidth = int(1000*sf)
textsizevalue = int(20*sf)
thumbnailsize = 1.5*sf

# We have a few rows of buttons at the bottom of the canvas, and their exact
# positions depend on the number of images. 
if (number_images <= 6):
	toprow_y = 720
	bottomrow_y = 790
if ((number_images > 6) & (number_images <= 12)):
	toprow_y = 870
	bottomrow_y = 940
if ((number_images > 12) & (number_images <= 18)):
	toprow_y = 1020
	bottomrow_y = 1060
if ((number_images > 18) & (number_images <= 32)):
	toprow_y = 1020
	bottomrow_y = 1060
	thumbnailsize = 1.2*sf

# These do not change depending on the number of images we have.
eazy_positionx, eazy_positiony = 500*sf, 230*sf
eazytext_positionx, eazytext_positiony = 350*sf, 70*sf
beagle_positionx, beagle_positiony = 1500*sf, 350*sf
beagletext_positionx, beagletext_positiony = 1110*sf, 70*sf#1300*sf, 70*sf

# Open up the photometric catalog
fitsinput = fits.open(input_photometry)
ID_values = fitsinput[1].data['ID'].astype('int')
RA_values = fitsinput[1].data['RA']
DEC_values = fitsinput[1].data['DEC']
number_objects = len(ID_values)

image_flux_value_cat = np.zeros([number_objects, number_image_filters])
image_flux_value_err_cat = np.zeros([number_objects, number_image_filters])
SNR_values = np.zeros([number_objects, number_image_filters])

for j in range(0, number_image_filters):
	image_flux_value_cat[:,j] = fitsinput[1].data[all_images_filter_name[j]]
	image_flux_value_err_cat[:,j] = fitsinput[1].data[all_images_filter_name[j]+'_err']
	SNR_values[:,j] = image_flux_value_cat[:,j] / image_flux_value_err_cat[:,j]

number_input_objects = len(ID_values)
ID_iterator = 0

# Decide whether or not the user requested an ID number or an id number list
if (args.id_number):
	ID_list = ID_values
	ID_list_indices = np.arange(len(ID_values), dtype = int)
	current_id = int(args.id_number)
	current_index = np.where(ID_values == current_id)[0][0]
	ID_iterator = current_index
	if (args.id_number_list):
		print("You can't specify an individual ID and a list, ignoring the list.")
	if (args.idarglist):
		print("You can't specify an individual ID and a list, ignoring the list.")
	
if not (args.id_number):
	if not (args.id_number_list):
		ID_list = ID_values
		ID_list_indices = np.arange(len(ID_values), dtype = int)
		current_index = ID_list_indices[ID_iterator]
		current_id = ID_list[current_index]

	if (args.id_number_list):
		ID_input_file = np.loadtxt(args.id_number_list)
		if (len(ID_input_file.shape) > 1):
			ID_numbers_to_view = ID_input_file[:,0].astype(int)
		else:
			ID_numbers_to_view = ID_input_file.astype(int)
		number_id_list = len(ID_numbers_to_view)
	
		# Set up index array for 
		ID_list_indices = np.zeros(number_id_list, dtype = int)
		for x in range(0, number_id_list):
			ID_list_indices[x] = np.where(ID_values == ID_numbers_to_view[x])[0]
	
		ID_list = ID_numbers_to_view
		current_index = ID_list_indices[ID_iterator]
		current_id = ID_values[current_index]

	if (args.idarglist):
		ID_numbers_to_view = np.array(ast.literal_eval(args.idarglist))
	
		number_id_list = len(ID_numbers_to_view)
	
		# Set up index array for 
		ID_list_indices = np.zeros(number_id_list, dtype = int)
		for x in range(0, number_id_list):
			ID_list_indices[x] = np.where(ID_values == ID_numbers_to_view[x])[0]
	
		ID_list = ID_numbers_to_view
		current_index = ID_list_indices[ID_iterator]
		current_id = ID_values[current_index]


# Create the notes array
notes_values = np.array([''], dtype = 'object')
for x in range(0, number_input_objects-1):
	notes_values = np.append(notes_values, ['']) 

# Create the flag arrays
highZflag_array = np.zeros(number_input_objects, dtype = 'int')
badfitflag_array = np.zeros(number_input_objects, dtype = 'int')
baddataflag_array = np.zeros(number_input_objects, dtype = 'int')

# So, now, there are three arrays:
#   ID_list (which is the list of the ID values that will be viewed)
#   ID_list_indices (which is the indices for the ID values that will be viewed)

# There is also a:
#   current_id - the current ID number being displayed
#   current_index - the current index in the full photometric array 
#   ID_iterator - an iterator that keeps track of the index of the ID_list_indices array
#                 is currently being shown



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Now, everything is set up, so let's start creating the GUI

# Start by creating the GUI as root
root=Tk()
root.wm_title("JADESView")

# Create the canvas 
canvas=Canvas(root, height=canvasheight, width=canvaswidth, bg="#ffffff")

# Plot the EAZY SED
image = Image.open(EAZY_files+str(current_id)+"_EAZY_SED.png")

# Crop out the thumbnails
image = cropEAZY(image)

photo = resizeimage(image)
item4 = canvas.create_image(eazy_positionx, eazy_positiony, image=photo)
Label(root, text="EAZY FIT", font=('helvetica', int(textsizevalue*1.5))).place(x=eazytext_positionx, y = eazytext_positiony)

# Plot the BEAGLE SED
new_image = Image.open(BEAGLE_files+str(current_id)+"_BEAGLE_SED.png")
new_photo = resizeimage(new_image)
item5 = canvas.create_image(beagle_positionx, beagle_positiony, image=new_photo)
Label(root, text="BEAGLE FIT", font=('helvetica', int(textsizevalue*1.5))).place(x=beagletext_positionx, y = beagletext_positiony)
	
canvas.pack(side = TOP, expand=True, fill=BOTH)

# Plot the thumbnails
fig_photo_objects = np.empty(0, dtype = 'object')
fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, current_id, current_index, defaultstretch)


# # # # # # # # # # # #
# Flag Object Buttons

# Create the Bad Fit Flag
btn1 = Button(root, text = 'Bad Fit', bd = '5', command = badfit)
btn1.config(height = int(2*sf), width = int(13*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue), padx = 3, pady = 3)
btn1.place(x = 600*sf, y = toprow_y*sf)

# Create the High Redshift Flag Button
btn1 = Button(root, text = 'High Redshift', bd = '5', command = highz)
btn1.config(height = int(2*sf), width = int(13*sf), fg='red', highlightbackground='white', font=('helvetica', textsizevalue), padx = 20, pady = 3)
btn1.place(x = 785*sf, y = toprow_y*sf)

# Create the Bad Data Flag
btn1 = Button(root, text = 'Bad Data', bd = '5', command = baddata)
btn1.config(height = int(2*sf), width = int(13*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue), padx = 3, pady = 3)
btn1.place(x = 1000*sf, y = toprow_y*sf)


# # # # # # # # # # # #
# Move to New Object Buttons


# Create the Previous Object Button
btn3 = Button(root, text = 'Previous Object', bd = '5', command = previousobject)  
btn3.config(height = int(2*sf), width = int(20*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue), padx = 3, pady = 3)
btn3.place(x = 600*sf, y = bottomrow_y*sf)


# Create the Next Object Button
btn2 = Button(root, text = 'Next Object', bd = '5', command = nextobject)
btn2.config(height = int(2*sf), width = int(20*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue), padx = 3, pady = 3)
btn2.place(x = 917*sf, y = bottomrow_y*sf)

if ((args.id_number_list is None) & (args.idarglist is None)):
	# Create the Object Entry Field and Button
	Label(root, text="Display Object: ", font=('helvetica', textsizevalue), fg="#000000", bg="#ffffff").place(x=1220*sf, y = (bottomrow_y+10.0)*sf)
	e1 = Entry(root, width = int(5*sf), font=('helvetica', textsizevalue), fg="#000000", bg="#ffffff")
	e1.place(x = 1370*sf, y = (bottomrow_y+6.0)*sf)

	btn9 = Button(root, text = 'Go', bd = '5', command = gotoobject)  
	btn9.config(height = int(1*sf), width = int(4*sf), fg='blue', highlightbackground='white', font=('helvetica', textsizevalue))
	#btn2.pack(side = 'bottom')
	btn9.place(x = 1470*sf, y = (bottomrow_y+9.0)*sf)


# # # # # # # # # # # #
# Quit Button

btn4 = Button(root, text = 'Quit', bd = '5', command = save_destroy)  
btn4.config(height = int(2*sf), width = int(10*sf), fg='red', highlightbackground='white', font=('helvetica', textsizevalue))
btn4.place(x = 1850*sf, y = (bottomrow_y+10.0)*sf)

# # # # # # # # # # # #
# Save Canvas Button

btn4 = Button(root, text = 'Save Canvas', bd = '5', command = save_canvas)  
btn4.config(height = int(2*sf), width = int(15*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue))
btn4.place(x = 1650*sf, y = (bottomrow_y+10.0)*sf)
 

# # # # # # # # # # # #
# Image Stretch Buttons

Label(root, text="Stretch", font = ('helvetica', int(20*sf)), fg="#000000", bg='#ffffff').place(x=20*sf, y = (bottomrow_y+ 10.0)*sf)

# Create the LinearStretch Button
btn5 = Button(root, text = 'Linear', bd = '5', command = linearstretch)  
if (defaultstretch == 'LinearStretch'):
	btn5.config(height = int(2*sf), width = int(10*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue))
else:
	btn5.config(height = int(2*sf), width = int(10*sf), fg='grey', highlightbackground='white', font=('helvetica', textsizevalue))
btn5.place(x = 100*sf, y = (bottomrow_y+10.0)*sf)


# Create the LogStretch Button
btn6 = Button(root, text = 'Log', bd = '5', command = logstretch)  
if (defaultstretch == 'LogStretch'):
	btn6.config(height = int(2*sf), width = int(10*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue))
else:
	btn6.config(height = int(2*sf), width = int(10*sf), fg='grey', highlightbackground='white', font=('helvetica', textsizevalue))
btn6.place(x = 250*sf, y = (bottomrow_y+10.0)*sf)


# Create the Asinh Button
btn7 = Button(root, text = 'Asinh', bd = '5', command = asinhstretch)  
if (defaultstretch == 'LinearStretch'):
	btn7.config(height = int(2*sf), width = int(10*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue))
else:
	btn7.config(height = int(2*sf), width = int(10*sf), fg='grey', highlightbackground='white', font=('helvetica', textsizevalue))
btn7.place(x = 400*sf, y = (bottomrow_y+10.0)*sf)


# Create the Notes Field
Label(root, text="Notes", font = "Helvetica 20", fg="#000000", bg="#ffffff").place(x=1220*sf, y = (toprow_y+5.0)*sf)
e2 = Entry(root, width = int(50*sf), font=('helvetica', textsizevalue), fg="#000000", bg="#ffffff")
e2.place(x = 1300*sf, y = toprow_y*sf)
e2.insert(0, notes_values[current_index])

# Create the RA and DEC size field 
Label(root, text="RA/DEC size", font=('helvetica', textsizevalue), fg="#000000", bg="#ffffff").place(x=20*sf, y = (toprow_y+15.0)*sf)
e3 = Entry(root, width = int(10*sf), font=('helvetica', textsizevalue), fg="#000000", bg="#ffffff")
e3.place(x = 150*sf, y = (toprow_y+10.0)*sf)
e3.insert(0, str(ra_dec_size_value))
Label(root, text="arcseconds", font=('helvetica', textsizevalue), fg="#000000", bg="#ffffff").place(x=280*sf, y = (toprow_y+15.0)*sf)
btn8 = Button(root, text = 'Change', bd = '5', command = changeradecsize)  
btn8.config(height = 1, width = int(10*sf), fg='blue', highlightbackground = 'white', font=('helvetica', textsizevalue))
btn8.place(x = 400*sf, y = (toprow_y+15.0)*sf)

root.mainloop()
