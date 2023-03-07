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

import webbrowser
def callback(url):
    webbrowser.open_new(url)

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

# Currently, for testing JADESView, we have z_spec
use_zspec = True

def getEAZYimage(ID):
	start_time = time.time()
	EAZY_file_name = EAZY_files+str(ID)+'_EAZY_SED.png'

	if (EAZY_file_name.startswith('http')):
		response = requests.get(EAZY_file_name, auth=HTTPBasicAuth(fenrir_username, fenrir_password))
		image = Image.open(BytesIO(response.content))
	else:
		image = Image.open(EAZY_file_name)
		
	end_time = time.time()
	if (timer_verbose):
		print("Fetching the EAZY image: " +str(end_time - start_time))

	return image

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

def highz_next():
	global current_index
	global ID_iterator
	global ID_list
	global highZflag_array
	
	highZflag_array[current_index] = 1
	current_id = ID_list[ID_iterator]
	print("Object "+str(current_id)+" is a high-redshift candidate. Moving on.")
	nextobject()

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

def baddata_next():
	global current_index
	global ID_iterator
	global ID_list
	global baddataflag_array
	
	baddataflag_array[current_index] = 1
	current_id = ID_list[ID_iterator]
	print("Object "+str(current_id)+" object has bad data. Moving on.")
	nextobject()


def update_eazy_text(current_id, eazy_results_IDs, eazy_results_zpeak):
	eazy_z_peak = getfile_value(current_id, eazy_results_IDs, eazy_results_zpeak, 4)
	eazy_z_a = getfile_value(current_id, eazy_results_IDs, eazy_results_za, 4)
	eazy_l68 = getfile_value(current_id, eazy_results_IDs, eazy_results_zl68, 4)
	eazy_u68 = getfile_value(current_id, eazy_results_IDs, eazy_results_zu68, 4)

	eazy_label_zpeak.configure(text="z_EAZY, peak = "+str(eazy_z_peak)+" ("+str(eazy_l68)+" - "+str(eazy_u68)+")")  
	eazy_label_za.configure(text="z_EAZY, a = "+str(eazy_z_a))  

def refresh_object():
	global e2
	global ID_iterator
	global current_index
	global ID_list
	global ID_list_indices
	global photo
	global new_photo
	global item4
	global canvas   
	global fig_photo_objects
	global defaultstretch

	global ID_values
	global RA_values
	global DEC_values

	global eazy_positionx, eazy_positiony
	global eazytext_positionx, eazytext_positiony
	
	notes_values[current_index] = e2.get()
	e2.delete(0,END)
	
	if (item4 is not None):
		canvas.delete(item4)
		
	current_index = ID_list_indices[ID_iterator]
	current_id = ID_list[ID_iterator]
	e2.insert(0, notes_values[current_index])

	#image = Image.open(EAZY_files+str(current_id)+"_EAZY_SED.png")
	if (EAZY_plots_exist == True):
		image = getEAZYimage(current_id)
		start_time = time.time()
		image = cropEAZY(image)
		end_time = time.time()
		if (timer_verbose):
			print("Cropping the EAZY image: " +str(end_time - start_time))
		start_time = time.time()
		photo = resizeimage(image)
		end_time = time.time()
		if (timer_verbose):
			print("Resizing the EAZY image: " +str(end_time - start_time))
		start_time = time.time()
		item4 = canvas.create_image(eazy_positionx, eazy_positiony, image=photo)
		end_time = time.time()
		if (timer_verbose):
			print("Creating the EAZY canvas: " +str(end_time - start_time))
	
	object_label.configure(text="Object "+str(current_id))  

	if (EAZY_results_file_exists):
		update_eazy_text(current_id, eazy_results_IDs, eazy_results_zpeak)

	if (fitsmap_link_exists == True):
		final_fitsmap_link = get_fitsmap_link(fitsmap_link, current_id, ID_values, RA_values, DEC_values)
		link.bind("<Button-1>", lambda e: callback(final_fitsmap_link))

	save_files()
		

def nextobject():
	global e2
	global ID_iterator
	global current_index
	global ID_list
	global ID_list_indices
	global photo
	global new_photo
	global item4
	global canvas   
	global fig_photo_objects
	global defaultstretch

	global ID_values
	global RA_values
	global DEC_values

	global eazy_positionx, eazy_positiony
	global eazytext_positionx, eazytext_positiony
	
	notes_values[current_index] = e2.get()
	e2.delete(0,END)

	if (ID_iterator < len(ID_list)-1):
		ID_iterator = ID_iterator+1
#	else:
#		len(ID_list)-1
#		save_destroy()
	
	if (item4 is not None):
		canvas.delete(item4)


		
	current_index = ID_list_indices[ID_iterator]
	current_id = ID_list[ID_iterator]
	e2.insert(0, notes_values[current_index])

	#image = Image.open(EAZY_files+str(current_id)+"_EAZY_SED.png")
	if (EAZY_plots_exist == True):
		image = getEAZYimage(current_id)
		start_time = time.time()
		image = cropEAZY(image)
		end_time = time.time()
		if (timer_verbose):
			print("Cropping the EAZY image: " +str(end_time - start_time))
		start_time = time.time()
		photo = resizeimage(image)
		end_time = time.time()
		if (timer_verbose):
			print("Resizing the EAZY image: " +str(end_time - start_time))
		start_time = time.time()
		item4 = canvas.create_image(eazy_positionx, eazy_positiony, image=photo)
		end_time = time.time()
		if (timer_verbose):
			print("Creating the EAZY canvas: " +str(end_time - start_time))
	
	object_label.configure(text="Object "+str(current_id))  

	if (EAZY_results_file_exists):
		update_eazy_text(current_id, eazy_results_IDs, eazy_results_zpeak)

	if (fitsmap_link_exists == True):
		final_fitsmap_link = get_fitsmap_link(fitsmap_link, current_id, ID_values, RA_values, DEC_values)
		link.bind("<Button-1>", lambda e: callback(final_fitsmap_link))

	save_files()


def previousobject():
	global ID_iterator
	global current_index
	global e2
	global ID_list
	global ID_list_indices
	global photo
	global new_photo
	global item4
	global canvas   
	global fig_photo_objects
	global defaultstretch

	global ID_values
	global RA_values
	global DEC_values

	global eazy_positionx, eazy_positiony
	global eazytext_positionx, eazytext_positiony

	notes_values[current_index] = e2.get()
	e2.delete(0,END)

	if (ID_iterator > 0):
		ID_iterator = ID_iterator-1
	else:
		ID_iterator = 0
	
	current_index = ID_list_indices[ID_iterator]
	current_id = ID_list[ID_iterator]
	e2.insert(0, notes_values[current_index])

	if (EAZY_plots_exist == True):
		canvas.delete(item4)
		#image = Image.open(EAZY_files+str(current_id)+"_EAZY_SED.png")
		image = getEAZYimage(current_id)
		image = cropEAZY(image)
		photo = resizeimage(image)
		item4 = canvas.create_image(eazy_positionx, eazy_positiony, image=photo)

	object_label.configure(text="Object "+str(current_id))  

	if (EAZY_results_file_exists):
		update_eazy_text(current_id, eazy_results_IDs, eazy_results_zpeak)

	if (fitsmap_link_exists == True):
		final_fitsmap_link = get_fitsmap_link(fitsmap_link, current_id, ID_values, RA_values, DEC_values)
		link.bind("<Button-1>", lambda e: callback(final_fitsmap_link))

	save_files()

def gotoobject():

	global ID_iterator
	global current_index
	global e2
	global ID_list
	global ID_list_indices
	global photo
	global new_photo
	global item4
	global canvas   
	global fig_photo_objects
	global defaultstretch

	global ID_values
	global RA_values
	global DEC_values

	global eazy_positionx, eazy_positiony
	global eazytext_positionx, eazytext_positiony

	notes_values[current_index] = e2.get()
	e2.delete(0,END)

	if (e1.get().isdigit() == True):

		#current_id = int(e1.get())
		ID_iterator = np.where(ID_list == int(e1.get()))[0][0]

		current_index = ID_list_indices[ID_iterator]
		current_id = ID_list[ID_iterator]
		e2.insert(0, notes_values[current_index])
	
		#image = Image.open(EAZY_files+str(current_id)+"_EAZY_SED.png")
		if (EAZY_plots_exist == True):
			canvas.delete(item4)
			image = getEAZYimage(current_id)
			image = cropEAZY(image)
			photo = resizeimage(image)
			item4 = canvas.create_image(eazy_positionx, eazy_positiony, image=photo)
		
		object_label.configure(text="Object "+str(current_id))  

		if (EAZY_results_file_exists):
			update_eazy_text(current_id, eazy_results_IDs, eazy_results_zpeak)

	else:
		print("That's not a valid ID number.")

	if (fitsmap_link_exists == True):
		final_fitsmap_link = get_fitsmap_link(fitsmap_link, current_id, ID_values, RA_values, DEC_values)
		link.bind("<Button-1>", lambda e: callback(final_fitsmap_link))

	save_files()


# This will remove the thumbnails, for future work
def cropEAZY(img):

	#output_image = img.crop((0, 0, 3300, 1600))
	#output_image = img.crop((0, 0, 3300, 1480))

	return img

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

def save_files():
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
	w.write('#ID  Flag   Notes \n')
	for z in range(0, len(ID_values)):
		output_flag = 4
#		if (notes_values[z] != ''):
		if (highZflag_array[z] == 1):
			output_flag = 1
		if (badfitflag_array[z] == 1):
			output_flag = 2
		if (baddataflag_array[z] == 1):
			output_flag = 3
		if (output_flag < 4):
			w.write(str(ID_values[z])+' '+str(output_flag)+' '+str(notes_values[z])+'\n')
	w.close()


def save_destroy():

	save_files()
	quit()
	#root.destroy()


def getfile_value(current_id, results_IDs, results_values, round_value):
	find_value_index = np.where(results_IDs == current_id)[0]
	if (len(find_value_index) > 0):
		return round(results_values[find_value_index[0]],round_value)
	else:
		return -9999

def getfile_true_or_false(current_id, results_IDs, results_values):
	find_value_index = np.where(results_IDs == current_id)[0]
	if (len(find_value_index) > 0):
		return results_values[find_value_index[0]] 
	else:
		return False

def get_fitsmap_link(fitsmap_link, objID, ID_values, RA_values, DEC_values):
	
	RA_object = RA_values[np.where(ID_values == objID)[0][0]]
	DEC_object = DEC_values[np.where(ID_values == objID)[0][0]]
	
	fitsmap_link = fitsmap_link.replace("RA_REPLACE", str(RA_object))
	fitsmap_link = fitsmap_link.replace("DEC_REPLACE", str(DEC_object))
	return fitsmap_link

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

# command line argument list of objects
parser.add_argument(
  '-idarglist',
  help="Command line argument list of objects",
  action="store",
  type=str,
  dest="idarglist",
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

timer_verbose = False
if (args.tverb):
	timer_verbose = True


# Right now, the default canvaswidth is 2000. 
canvaswidth = 2000

# I have to set these as false unless the file is specified in the input file
EAZY_plots_exist = False
EAZY_results_file_exists = False

fitsmap_link_exists = False
# Read in the various input values from the input file. 
input_lines = np.loadtxt(JADESView_input_file, dtype='str')
number_input_lines = len(input_lines[:,0])
for i in range(0, number_input_lines):
	if (input_lines[i,0] == 'input_photometry'):
		input_photometry = input_lines[i,1]
	if (input_lines[i,0] == 'EAZY_files'):
		EAZY_files = input_lines[i,1]
		EAZY_plots_exist = True
	if (input_lines[i,0] == 'EAZY_results'):
		EAZY_results_file = input_lines[i,1]
		EAZY_results_file_exists = True
		EAZY_results_file_exists = os.path.exists(EAZY_results_file)
	if (input_lines[i,0] == 'output_flags_file'):
		output_flags_file = input_lines[i,1]
	if (input_lines[i,0] == 'output_notes_file'):
		output_notes_file = input_lines[i,1]
	if (input_lines[i,0] == 'canvaswidth'):
		canvaswidth = float(input_lines[i,1])
	if (input_lines[i,0] == 'fenrir_username'):
		fenrir_username = input_lines[i,1]
	if (input_lines[i,0] == 'fenrir_password'):
		fenrir_password = input_lines[i,1]
	if (input_lines[i,0] == 'fitsmap_link'):
		fitsmap_link = input_lines[i,1]
		fitsmap_link_exists = True

#base64string = base64.b64encode('%s:%s' % (fenrir_username, fenrir_password))

# # # # # # # # # # # # # # # # # # 
# Let's open up all the input files



sf = canvaswidth / 2000.0 # This is the "shrinkfactor" by which all of the canvas
                          # element positions and sizes are shrunk or expanded. I 
                          # RECOGNIZE THAT I SHOULD PUT THINGS ON A GRID, BUT THAT
                          # WILL COME IN A FUTURE UPDATE, OK

# The fontsize depends on this
fontsize = str(int(20*sf))

canvasheight = (canvaswidth*(1.0 / 1.3))  # For showing the various results on the figure, 
										   # I lock everything to a 1.8:1 aspect ratio

baseplotwidth = int(1900*sf)
textsizevalue = int(20*sf)
thumbnailsize = 1.5*sf

toprow_y = 1020
bottomrow_y = 1060

# These do not change depending on the number of images we have.
eazy_positionx, eazy_positiony = 1000*sf, 800*sf#245*sf
eazytext_positionx, eazytext_positiony = 820*sf, 10*sf
objectID_positionx, objectID_positiony = 20*sf, 10*sf#1300*sf, 70*sf

# Open up the photometric catalog
fitsinput = fits.open(input_photometry)
ID_values = fitsinput[1].data['ID'].astype('int')
RA_values = fitsinput[1].data['RA']
DEC_values = fitsinput[1].data['DEC']
number_objects = len(ID_values)

number_input_objects = len(ID_values)
ID_iterator = 0

if (EAZY_results_file_exists):
	if (EAZY_results_file.startswith('http')):
		response = requests.get(EAZY_results_file, auth=HTTPBasicAuth(fenrir_username, fenrir_password))
		eazy_fits_file = BytesIO(response.content)
	else:
		eazy_fits_file = EAZY_results_file
	eazy_results_fits = fits.open(eazy_fits_file)
	eazy_results_IDs = eazy_results_fits[1].data['ID'].astype('int')
	eazy_results_zpeak = eazy_results_fits[1].data['z_peak']
	eazy_results_za = eazy_results_fits[1].data['z_a']
	eazy_results_zl68 = eazy_results_fits[1].data['l68']
	eazy_results_zu68 = eazy_results_fits[1].data['u68']
	#eazy_results_zl95 = eazy_results_fits[1].data['l95']
	#eazy_results_zu95 = eazy_results_fits[1].data['u95']

	
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
			#print(np.where(ID_values == ID_numbers_to_view[x])[0])
			#print(len(np.where(ID_values == ID_numbers_to_view[x])[0]))
			if (len(np.where(ID_values == ID_numbers_to_view[x])[0]) > 0):
				ID_list_indices[x] = np.where(ID_values == ID_numbers_to_view[x])[0]
			else:
				sys.exit("Object "+str(ID_numbers_to_view[x])+" does not appear in this catalog. Exiting.")
					
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
#image = Image.open(EAZY_files+str(current_id)+"_EAZY_SED.png")
image = getEAZYimage(current_id)

# Put the object label 
object_label = Label(root, text="Object "+str(current_id), font = "Helvetica "+str(int(textsizevalue*1.5)), fg="black", bg="white")
object_label.place(x=objectID_positionx, y = objectID_positiony)

# Crop out the thumbnails
if (EAZY_plots_exist == True):
	image = cropEAZY(image)
	photo = resizeimage(image)
	item4 = canvas.create_image(eazy_positionx, eazy_positiony, image=photo)
	#Label(root, text="EAZY FIT", fg='black', bg='white', font=('helvetica', int(textsizevalue*1.5))).place(x=eazytext_positionx, y = eazytext_positiony)
else:
	item4 = None

canvas.pack(side = TOP, expand=True, fill=BOTH)


# # # # # # # # # # # # # # 
# Place Labels with Redshift 

# #711c91, #ea00d9, #0abdc6, #133e7c, #091833

# A delineation line. 
#redshift_separator = canvas.create_rectangle(1100*sf, (toprow_y-320.0)*sf, 1940*sf, (toprow_y-310.0)*sf, outline="#0abdc6", fill="#0abdc6", tags="separator")


# Make the EAZY redshift label
if (EAZY_results_file_exists == True):
	eazy_z_peak = getfile_value(current_id, eazy_results_IDs, eazy_results_zpeak, 4)
	eazy_z_a = getfile_value(current_id, eazy_results_IDs, eazy_results_za, 4)
	eazy_l68 = getfile_value(current_id, eazy_results_IDs, eazy_results_zl68, 4)
	eazy_u68 = getfile_value(current_id, eazy_results_IDs, eazy_results_zu68, 4)

	eazy_label_za = Label(root, text="z_EAZY, a = "+str(eazy_z_a)+" ("+str(eazy_l68)+" - "+str(eazy_u68)+")", font = "Helvetica "+str(textsizevalue), fg="#133e7c", bg="#ffffff")
	eazy_label_za.place(x=1100*sf, y = 1250*sf)
	eazy_label_zpeak = Label(root, text="z_EAZY, peak = "+str(eazy_z_peak), font = "Helvetica "+str(textsizevalue), fg="#133e7c", bg="#ffffff")
	eazy_label_zpeak.place(x=1100*sf, y = 1300*sf)



# # # # # # # # # # # #
# Flag Object Buttons

# Create the Bad Fit Flag
btn1 = Button(root, text = 'Bad Fit', bd = '5', command = badfit)
btn1.config(height = int(2*sf), width = int(13*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue), padx = 3, pady = 3)
btn1.place(x = 600*sf, y = 1400*sf)

# Create the High Redshift Flag Button
btn1 = Button(root, text = 'High Redshift', bd = '5', command = highz)
btn1.config(height = int(2*sf), width = int(11*sf), fg='red', highlightbackground='white', font=('helvetica', textsizevalue), padx = 20, pady = 3)
btn1.place(x = 765*sf, y = 1400*sf)

# Create the High Redshift, Next Object Flag Button
btn1 = Button(root, text = 'High Redshift, Next Object', bd = '5', command = highz_next)
btn1.config(height = int(2*sf), width = int(20*sf), fg='red', highlightbackground='white', font=('helvetica', textsizevalue), padx = 20, pady = 3)
btn1.place(x = 1300*sf, y = 1380*sf)

# Create the Bad Data Flag
btn1 = Button(root, text = 'Bad Data', bd = '5', command = baddata)
btn1.config(height = int(2*sf), width = int(13*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue), padx = 3, pady = 3)
btn1.place(x = 975*sf, y = 1400*sf)

# Create the Bad Data Flag, Next Object
btn1 = Button(root, text = 'Bad Data, Next Object', bd = '5', command = baddata_next)
btn1.config(height = int(2*sf), width = int(20*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue), padx = 3, pady = 3)
btn1.place(x = 1600*sf, y = 1380*sf)

# refresh_object
btn1 = Button(root, text = 'Refresh', bd = '5', command = refresh_object)
btn1.config(height = int(2*sf), width = int(20*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue), padx = 3, pady = 3)
btn1.place(x = 300*sf, y = 1425*sf)

# # # # # # # # # # # #
# Move to New Object Buttons


# Create the Previous Object Button
btn3 = Button(root, text = 'Previous Object', bd = '5', command = previousobject)  
btn3.config(height = int(2*sf), width = int(20*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue), padx = 3, pady = 3)
btn3.place(x = 600*sf, y = 1450*sf)


# Create the Next Object Button
btn2 = Button(root, text = 'Next Object', bd = '5', command = nextobject)
btn2.config(height = int(2*sf), width = int(20*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue), padx = 3, pady = 3)
btn2.place(x = 917*sf, y = 1450*sf)

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
btn4.place(x = 1850*sf, y = 1450*sf)

if (fitsmap_link_exists == True):
	link = Label(root, text="Fitsmap Link",font=('Helveticabold', 15), fg="blue", cursor="hand2")
	link.pack()
	final_fitsmap_link = get_fitsmap_link(fitsmap_link, current_id, ID_values, RA_values, DEC_values)
	link.bind("<Button-1>", lambda e: callback(final_fitsmap_link))


# # # # # # # # # # # #
# Save Canvas Button

# save_canvas_imagegrab()
#btn4 = Button(root, text = 'Save Canvas', bd = '5', command = save_canvas)  
#btn4.config(height = int(2*sf), width = int(15*sf), fg='black', highlightbackground='white', font=('helvetica', textsizevalue))
#btn4.place(x = 1645*sf, y = (bottomrow_y+10.0)*sf)


# # # # # # # # #  
# Notes / RA/DEC

# Create the Notes Field
Label(root, text="Notes", font = "Helvetica 20", fg="#000000", bg="#ffffff").place(x=1200*sf, y = 1450*sf)
e2 = Entry(root, width = int(70*sf), font=('helvetica', textsizevalue), fg="#000000", bg="#ffffff")
e2.place(x = 1300*sf, y = 1450*sf)
e2.insert(0, notes_values[current_index])


# # # # # # # # #  
# Alternate Fits 


root.mainloop()