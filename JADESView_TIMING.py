import os
import ast
import sys
import math
import time
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

# Right now, everything is kind of built with this as the size of things
canvasheight = 1000
canvaswidth = 2000

eazy_positionx, eazy_positiony = 500, 230
eazytext_positionx, eazytext_positiony = 350, 70
beagle_positionx, beagle_positiony = 1500, 350
beagletext_positionx, beagletext_positiony = 1300, 70

# The default stretch on the various images
defaultstretch = 'LogStretch'

# The default size of the various images
ra_dec_size_value = 2.0


def resizeimage(image):
	basewidth = 1000
	wpercent = (basewidth / float(image.size[0]))
	hsize = int((float(image.size[1]) * float(wpercent)))
	image = image.resize((basewidth, hsize), PIL.Image.ANTIALIAS)
	photo = ImageTk.PhotoImage(image)
	return photo
	

def highz():
	global current_index
	global ID_iterator
	global ID_list
	global highZflag_array
	
	highZflag_array[current_index] = 1
	current_id = ID_list[ID_iterator]
	print "Object "+str(current_id)+" is a high-redshift candidate."

def badfit():
	global current_index
	global badfitflag_array
	
	badfitflag_array[current_index] = 1
	
	print "Object "+str(current_id)+" has a bad fit."

def baddata():
	global current_index
	global baddataflag_array
	
	baddataflag_array[current_index] = 1
	print "Object "+str(current_id)+" object has bad data."


def nextobject():
	#print "Moving to the next object."
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
	else:
		ID_iterator = len(ID_list)-1
	
	current_index = ID_list_indices[ID_iterator]
	current_id = ID_list[ID_iterator]
	e2.insert(0, notes_values[current_index])

	#print id
	image = Image.open(EAZY_files+str(current_id)+"_EAZY_SED.png")
	image = cropEAZY(image)
	photo = resizeimage(image)
	item4 = canvas.create_image(eazy_positionx, eazy_positiony, image=photo)
	
	new_image = Image.open(BEAGLE_files+str(current_id)+"_BEAGLE_SED.png")
	new_photo = resizeimage(new_image)
	item5 = canvas.create_image(beagle_positionx, beagle_positiony, image=new_photo)

	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, current_id, current_index, defaultstretch)

def previousobject():
	#print "Moving to the next object."
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

	#print id
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
	
		#print id
		image = Image.open(EAZY_files+str(current_id)+"_EAZY_SED.png")
		image = cropEAZY(image)
		photo = resizeimage(image)
		item4 = canvas.create_image(eazy_positionx, eazy_positiony, image=photo)
		
		new_image = Image.open(BEAGLE_files+str(current_id)+"_BEAGLE_SED.png")
		new_photo = resizeimage(new_image)
		item5 = canvas.create_image(beagle_positionx, beagle_positiony, image=new_photo)
	
		fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, current_id, current_index, defaultstretch)
	else:
		print "That's not a valid ID number."


# This will remove the thumbnails, for future work
def cropEAZY(img):

	output_image = img.crop((0, 0, 3300, 1600))

	return output_image

def linearstretch():
	global ID_iterator
	global ID_list
	global ID_list_indices
	global canvas   
	global fig_photo_objects
	global defaultstretch
	global btn5
	global btn6
	global btn7

	btn5.config(height = 2, width = 10, fg='black', font=('helvetica', 20))
	btn6.config(height = 2, width = 10, fg='grey', font=('helvetica', 20))
	btn7.config(height = 2, width = 10, fg='grey', font=('helvetica', 20))

	defaultstretch = 'LinearStretch'	
	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, ID_list[ID_iterator], ID_list_indices[ID_iterator], defaultstretch)

def logstretch():
	global ID_iterator
	global ID_list
	global ID_list_indices
	global canvas   
	global fig_photo_objects
	global defaultstretch
	global btn5
	global btn6
	global btn7

	btn5.config(height = 2, width = 10, fg='grey', font=('helvetica', 20))
	btn6.config(height = 2, width = 10, fg='black', font=('helvetica', 20))
	btn7.config(height = 2, width = 10, fg='grey', font=('helvetica', 20))

	defaultstretch = 'LogStretch'
	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, ID_list[ID_iterator], ID_list_indices[ID_iterator], defaultstretch)

def asinhstretch():
	global ID_iterator
	global ID_list
	global ID_list_indices
	global canvas   
	global fig_photo_objects
	global defaultstretch
	global btn5
	global btn6
	global btn7

	btn5.config(height = 2, width = 10, fg='grey', font=('helvetica', 20))
	btn6.config(height = 2, width = 10, fg='grey', font=('helvetica', 20))
	btn7.config(height = 2, width = 10, fg='black', font=('helvetica', 20))

	defaultstretch = 'AsinhStretch'
	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, ID_list[ID_iterator], ID_list_indices[ID_iterator], defaultstretch)

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
    new_start_time = time.time()
    figure_canvas_agg = FigureCanvasAgg(figure)
    figure_canvas_agg.draw()
    figure_x, figure_y, figure_w, figure_h = figure.bbox.bounds
    figure_w, figure_h = int(figure_w), int(figure_h)
    photo = PhotoImage(master=canvas, width=figure_w, height=figure_h)
    print "      Running PhotoImage took %s seconds" % (time.time() - new_start_time)

    new_start_time = time.time()
    # Position: convert from top-left anchor to center anchor
    canvas.create_image(loc[0] + figure_w/2, loc[1] + figure_h/2, image=photo)
    print "       Running Canvas Create Image took %s seconds" % (time.time() - new_start_time)

    new_start_time = time.time()
    # Unfortunately, there's no accessor for the pointer to the native renderer
    tkagg.blit(photo, figure_canvas_agg.get_renderer()._renderer, colormode=2)
    print "        Running tkagg.blit took %s seconds" % (time.time() - new_start_time)

    # Return a handle which contains a reference to the photo object
    # which must be kept live or else the picture disappears
    return photo

def create_thumbnails(canvas, fig_photo_objects, id_value, id_value_index, stretch):
	#global ID_values
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
	start_time_full = time.time()
	for i in range(0, number_images):
		
		image = image_all[i].data
		image_hdu = image_hdu_all[i]
		image_wcs = image_wcs_all[i]
				
		if (all_images_filter_name[i] == 'HST_F814W'):
			image_wcs.sip = None
		
		if (image_flux_value_err_cat[idx_cat, i] > -9999):
			# Make the cutout
			
			start_time = time.time()
			image_cutout = Cutout2D(image, position, size, wcs=image_wcs)
			print "   Running Cutout2D took %s seconds" % (time.time() - start_time)

			# Create the wcs axes
			plt.clf()
			start_time = time.time()
			fig = plt.figure(figsize=(1.5,1.5))
			ax3 = fig.add_axes([0, 0, 1, 1], projection=image_cutout.wcs)
			ax3.text(0.51, 0.96, all_images_filter_name[i].split('_')[1], transform=ax3.transAxes, fontsize=12, fontweight='bold', ha='center', va='top', color = 'black')
			ax3.text(0.5, 0.95, all_images_filter_name[i].split('_')[1], transform=ax3.transAxes, fontsize=12, fontweight='bold', ha='center', va='top', color = 'white')
			if (SNR_values[idx_cat, i] > -100):
				ax3.text(0.96, 0.06, 'SNR = '+str(round(SNR_values[idx_cat, i],2)), transform=ax3.transAxes, fontsize=12, fontweight='bold', horizontalalignment='right', color = 'black')
				ax3.text(0.95, 0.05, 'SNR = '+str(round(SNR_values[idx_cat, i],2)), transform=ax3.transAxes, fontsize=12, fontweight='bold', horizontalalignment='right', color = 'white')
			else:
				ax3.text(0.96, 0.06, 'SNR < -100', transform=ax3.transAxes, fontsize=12, fontweight='bold', horizontalalignment='right', color = 'black')
				ax3.text(0.95, 0.05, 'SNR < -100', transform=ax3.transAxes, fontsize=12, fontweight='bold', horizontalalignment='right', color = 'white')
			print "    Plotting the text and SNR values took %s seconds" % (time.time() - start_time)
			
			# Set the color map
			plt.set_cmap('gray')
			
			start_time = time.time()
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
			print "     ImageNormalize took %s seconds" % (time.time() - start_time)
			
			start_time = time.time()
			if (indexerror == 0):
				ax3.imshow(thumbnail, origin = 'lower', aspect='equal', norm = norm)
			else:
				ax3.imshow(thumbnail, origin = 'lower', aspect='equal')
			print "      Running imshow took %s seconds" % (time.time() - start_time)
			
			start_time = time.time()
			if (i <= 5):
				fig_x, fig_y = 20+(175*i), 500
			if ((i > 5) & (i <= 11)):
				fig_x, fig_y = 20+(175*(i-6)), 675
			if ((i > 11) & (i <= 17)):
				fig_x, fig_y = 20+(175*(i-12)), 900
				
			# Keep this handle alive, or else figure will disappear
			fig_photo_objects = np.append(fig_photo_objects, draw_figure(canvas, fig, loc=(fig_x, fig_y)))
			plt.close('all')
			print "         Drawing the thumnails to the figures took %s seconds" % (time.time() - start_time)
	print "**** FULLY PLOTTING THUMBNAILS TOOK %s seconds" % (time.time() - start_time_full)

	return fig_photo_objects



def save_destroy():
	global ID_values
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

	notes_values[ID_iterator] = e2.get()

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

start_time = time.time()

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
	if (input_lines[i,0] == 'canvasheight'):
		canvasheight = float(input_lines[i,1])
	if (input_lines[i,0] == 'canvaswidth'):
		canvaswidth = float(input_lines[i,1])
	if (input_lines[i,0] == 'defaultstretch'):
		defaultstretch = input_lines[i,1]
	if (input_lines[i,0] == 'ra_dec_size_value'):
		ra_dec_size_value = float(input_lines[i,1])
	
print "Opening up the input file took %s seconds" % (time.time() - start_time)


# # # # # # # # # # # # # # # # # # 
# Let's open up all the input files


start_time = time.time()
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
	print "  Opening up image: "+all_image_paths[i]
	image_all = np.append(image_all, fits.open(all_image_paths[i]))
	image_hdu_all = np.append(image_hdu_all, fits.open(all_image_paths[i])[0])
	image_wcs_all = np.append(image_wcs_all, WCS(image_hdu_all[i].header))

print "Opening up the individual image files took %s seconds" % (time.time() - start_time)

start_time = time.time()

# Open up the photometric catalog
fitsinput = fits.open(input_photometry)
ID_values = fitsinput[1].data['ID']
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

print "Opening up the photometric catalog took %s seconds" % (time.time() - start_time)

start_time = time.time()

number_input_objects = len(ID_values)
ID_iterator = 0

# Decide whether or not the user requested an ID number or an id number list
if (args.id_number):
	ID_list = ID_values
	ID_list_indices = ID_values - 1
	current_id = int(args.id_number)
	current_index = np.where(ID_values == current_id)[0][0]
	ID_iterator = current_index
	if (args.id_number_list):
		print "You can't specify an individual ID and a list, ignoring the list."
	if (args.idarglist):
		print "You can't specify an individual ID and a list, ignoring the list."
	
if not (args.id_number):
	if not (args.id_number_list):
		ID_list = ID_values
		ID_list_indices = ID_values - 1
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

print "Setting up the ID list took %s seconds" % (time.time() - start_time)


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

start_time = time.time()

# Start by creating the GUI as root
root=Tk()
root.wm_title("JADESView")

# Create the canvas 
canvas=Canvas(root, height=canvasheight, width=canvaswidth)

print "Setting up the GUI took %s seconds" % (time.time() - start_time)

start_time = time.time()

# Plot the EAZY SED
image = Image.open(EAZY_files+str(current_id)+"_EAZY_SED.png")

# Crop out the thumbnails
image = cropEAZY(image)

photo = resizeimage(image)
item4 = canvas.create_image(eazy_positionx, eazy_positiony, image=photo)
Label(root, text="EAZY FIT", font = "Helvetica 30").place(x=eazytext_positionx, y = eazytext_positiony)

print "Plotting the initial EAZY file took %s seconds" % (time.time() - start_time)

start_time = time.time()

# Plot the BEAGLE SED
new_image = Image.open(BEAGLE_files+str(current_id)+"_BEAGLE_SED.png")
new_photo = resizeimage(new_image)
item5 = canvas.create_image(beagle_positionx, beagle_positiony, image=new_photo)
Label(root, text="BEAGLE FIT", font = "Helvetica 30").place(x=beagletext_positionx, y = beagletext_positiony)

print "Plotting the initial BEAGLE file took %s seconds" % (time.time() - start_time)
	

canvas.pack(side = TOP, expand=True, fill=BOTH)


start_time = time.time()

# Plot the thumbnails
fig_photo_objects = np.empty(0, dtype = 'object')
fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, current_id, current_index, defaultstretch)

print "Creating the initial thumbnails took %s seconds" % (time.time() - start_time)

# # # # # # # # # # # #
# Flag Object Buttons

# Create the Bad Fit Flag
btn1 = Button(root, text = 'Bad Fit', bd = '5', command = badfit)
btn1.config(height = 2, width = 13, fg='black', font=('helvetica', 20))
btn1.place(x = 600, y = 870)

# Create the High Redshift Flag Button
btn1 = Button(root, text = 'High Redshift', bd = '5', command = highz)
btn1.config(height = 2, width = 13, fg='red', font=('helvetica', 20))
btn1.place(x = 800, y = 870)

# Create the Bad Data Flag
btn1 = Button(root, text = 'Bad Data', bd = '5', command = baddata)
btn1.config(height = 2, width = 13, fg='black', font=('helvetica', 20))
btn1.place(x = 1000, y = 870)


# # # # # # # # # # # #
# Move to New Object Buttons


# Create the Previous Object Button
btn3 = Button(root, text = 'Previous Object', bd = '5', command = previousobject)  
btn3.config(height = 2, width = 20, fg='black', font=('helvetica', 20))
btn3.place(x = 600, y = 940)

# Create the Next Object Button
btn2 = Button(root, text = 'Next Object', bd = '5', command = nextobject)
btn2.config(height = 2, width = 20, fg='black', font=('helvetica', 20))
btn2.place(x = 917, y = 940)

if ((args.id_number_list is None) & (args.idarglist is None)):
	# Create the Object Entry Field and Button
	Label(root, text="Display Object: ", font = "Helvetica 20").place(x=1220, y = 950)
	e1 = Entry(root, width = 5, font = "Helvetica 20")
	e1.place(x = 1370, y = 946)
	
	btn9 = Button(root, text = 'Go', bd = '5', command = gotoobject)  
	btn9.config(height = 1, width = 4, fg='blue', font=('helvetica', 20))
	#btn2.pack(side = 'bottom')
	btn9.place(x = 1470, y = 949)


# # # # # # # # # # # #
# Quit Buttons

# Create the Quit Button
btn4 = Button(root, text = 'Quit', bd = '5', command = save_destroy)  
btn4.config(height = 2, width = 20, fg='grey', font=('helvetica', 20))
btn4.place(x = 1700, y = 940)

# # # # # # # # # # # #
# Image Stretch Buttons

Label(root, text="Stretch", font = "Helvetica 20").place(x=20, y = 950)

# Create the LinearStretch Button
btn5 = Button(root, text = 'Linear', bd = '5', command = linearstretch)  
if (defaultstretch == 'LinearStretch'):
	btn5.config(height = 2, width = 10, fg='black', font=('helvetica', 20))
else:
	btn5.config(height = 2, width = 10, fg='grey', font=('helvetica', 20))
btn5.place(x = 100, y = 940)

# Create the LogStretch Button
btn6 = Button(root, text = 'Log', bd = '5', command = logstretch)  
if (defaultstretch == 'LogStretch'):
	btn6.config(height = 2, width = 10, fg='black', font=('helvetica', 20))
else:
	btn6.config(height = 2, width = 10, fg='grey', font=('helvetica', 20))
btn6.place(x = 250, y = 940)

# Create the Asinh Button
btn7 = Button(root, text = 'Asinh', bd = '5', command = asinhstretch)  
if (defaultstretch == 'AsinhStretch'):
	btn7.config(height = 2, width = 10, fg='black', font=('helvetica', 20))
else:
	btn7.config(height = 2, width = 10, fg='grey', font=('helvetica', 20))
btn7.place(x = 400, y = 940)




# Create the Notes Field
Label(root, text="Notes", font = "Helvetica 20").place(x=1220, y = 875)
e2 = Entry(root, width = 50, font = "Helvetica 20")
e2.place(x = 1300, y = 870)
e2.insert(0, notes_values[current_index])

# Create the RA and DEC size field 
Label(root, text="RA/DEC size", font = "Helvetica 20").place(x=20, y = 885)
e3 = Entry(root, width = 10, font = "Helvetica 20")
e3.place(x = 150, y = 880)
e3.insert(0, str(ra_dec_size_value))
Label(root, text="arcseconds", font = "Helvetica 20").place(x=280, y = 885)
btn8 = Button(root, text = 'Change', bd = '5', command = changeradecsize)  
btn8.config(height = 1, width = 10, fg='blue', font=('helvetica', 20))
btn8.place(x = 400, y = 885)



root.mainloop()
