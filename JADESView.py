import os
import sys
import math
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table

from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize, AsinhStretch, SinhStretch, LinearStretch)

import matplotlib.backends.tkagg as tkagg
from matplotlib.backends.backend_agg import FigureCanvasAgg

from Tkinter import *
import PIL
from PIL import ImageTk, Image

#from JADESView_functions import *

# Eventually this, and perhaps the arguments, will be parsed in a separate file. 
EAZY_files = '/Volumes/KNH_EXTERNAL/December_2019_DataChallenge/EAZY_run_1_13_20/F200W_plots_kron_f80/'
BEAGLE_files = '/Volumes/KNH_EXTERNAL/December_2019_DataChallenge/BEAGLE_output_1_13_20/F200W_kron_f80/'
all_filters_file_name = '/Users/knh/Desktop/NIRCam/photometric_redshifts/PhotoZReady/JADES_All_Filters.dat'
output_flags_file = 'Object_Flags.fits'
output_notes_file = 'Object_Notes.txt'

# Right now, everything is kind of built with this as the size of things
canvasheight = 1000
canvaswidth = 2000

eazy_positionx, eazy_positiony = 500, 230
eazytext_positionx, eazytext_positiony = 350, 70
beagle_positionx, beagle_positiony = 1500, 350
beagletext_positionx, beagletext_positiony = 1300, 70

# The default stretch on the various images
defaultstretch = 'AsinhStretch'

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
	global ID_list
	global highZflag_array
	
	highZflag_array[current_index] = 1
	current_id = ID_list[current_index]
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
	if (ID_iterator >= 0):
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

	global photo
	global new_photo
	global item4
	global item5
	global canvas   

	if (e1.get().isdigit() == True):
		canvas.delete(item4)
		canvas.delete(item5)
		id_value = 1
		id_value = int(e1.get())
		
		#print id
		image = Image.open(EAZY_files+str(id_value)+"_EAZY_SED.png")
		photo = resizeimage(image)
		item4 = canvas.create_image(500, 500, image=photo)
		
		new_image = Image.open(BEAGLE_files+str(id_value)+"_BEAGLE_SED.png")
		new_photo = resizeimage(new_image)
		item5 = canvas.create_image(1500, 500, image=new_photo)
	else:
		print "That's not a number."

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

	defaultstretch = 'LinearStretch'
	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, ID_list[ID_list_indices[ID_iterator]], ID_list_indices[ID_iterator], defaultstretch)

def sqrtstretch():
	global ID_iterator
	global ID_list
	global ID_list_indices
	global canvas   
	global fig_photo_objects
	global defaultstretch

	defaultstretch = 'SqrtStretch'
	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, ID_list[ID_list_indices[ID_iterator]], ID_list_indices[ID_iterator], defaultstretch)

def asinhstretch():
	global ID_iterator
	global ID_list
	global ID_list_indices
	global canvas   
	global fig_photo_objects
	global defaultstretch

	defaultstretch = 'AsinhStretch'
	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, ID_list[ID_list_indices[ID_iterator]], ID_list_indices[ID_iterator], defaultstretch)

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
	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, ID_list[ID_list_indices[ID_iterator]], ID_list_indices[ID_iterator], defaultstretch)


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
	
	for i in range(0, number_images):
		image = image_all[i].data
		image_hdu = image_hdu_all[i]
		image_wcs = image_wcs_all[i]
				
		if (all_images_filter_name[i] == 'HST_F814W'):
			image_wcs.sip = None
				
		if (image_flux_value_err_cat[idx_cat, i] > 0):
			# Make the cutout
			image_cutout = Cutout2D(image, position, size, wcs=image_wcs)
			
			# Create the wcs axes
			plt.clf()
			fig = plt.figure(figsize=(1.5,1.5))
			ax3 = fig.add_axes([0, 0, 1, 1], projection=image_cutout.wcs)
			ax3.text(0.5, 0.95, all_images_filter_name[i].split('_')[1], transform=ax3.transAxes, fontsize=12, fontweight='bold', ha='center', va='top', color = 'white')
			ax3.text(0.95, 0.05, 'SNR = '+str(round(SNR_values[idx_cat, i],2)), transform=ax3.transAxes, fontsize=12, fontweight='bold', horizontalalignment='right', color = 'white')
			# Set the color map
			plt.set_cmap('gray')
					
			# Normalize the image using the min-max interval and a square root stretch
			thumbnail = image_cutout.data
			if (stretch == 'AsinhStretch'):
				norm = ImageNormalize(thumbnail, interval=MinMaxInterval(), stretch=AsinhStretch())
			if (stretch == 'SqrtStretch'):
				norm = ImageNormalize(thumbnail, interval=MinMaxInterval(), stretch=SqrtStretch())
			if (stretch == 'LinearStretch'):
				norm = ImageNormalize(thumbnail, interval=MinMaxInterval(), stretch=LinearStretch())
			ax3.imshow(thumbnail, origin = 'lower', aspect='equal', norm = norm)
							
			if (i <= 5):
				fig_x, fig_y = 20+(175*i), 500
				#Label(root, text=all_images_filter_name[i].split('_')[1], font = "Helvetica 20").place(x=fig_x, y = fig_y)
			if ((i > 5) & (i <= 11)):
				fig_x, fig_y = 20+(175*(i-6)), 675
				#Label(root, text=all_images_filter_name[i].split('_')[1], font = "Helvetica 20").place(x=fig_x, y = fig_y)
			if ((i > 11) & (i <= 17)):
				fig_x, fig_y = 20+(175*(i-12)), 900
				#Label(root, text=all_images_filter_name[i].split('_')[1], font = "Helvetica 20").place(x=fig_x, y = fig_y)
				
			# Keep this handle alive, or else figure will disappear
			fig_photo_objects = np.append(fig_photo_objects, draw_figure(canvas, fig, loc=(fig_x, fig_y)))
			plt.close('all')

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

	#print notes_values
	quit()
	#root.destroy()


######################
# Required Arguments #
######################

parser = argparse.ArgumentParser()

# Input Photometry
parser.add_argument(
  '-input','--input_photometry',
  help="Input Photometry for analysis",
  action="store",
  type=str,
  dest="input_photometry",
  required=True
)

# Filters file
parser.add_argument(
  '-filters','--filters',
  help="Input file with filters that were used",
  action="store",
  type=str,
  dest="filters",
  required=True
)

# Input image file
parser.add_argument(
  '-ilist','--image_list',
  help="List of input images thumbnails are desired for",
  action="store",
  type=str,
  dest="image_list",
  required=False
)

######################
# Optional Arguments #
######################

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

args=parser.parse_args()

# # # # # # # # # # # # # # # # # # 
# Let's open up all the input files

# Open up the image list file
all_images_file_name = args.image_list
images_all_txt = np.loadtxt(all_images_file_name, dtype='str')
all_images_filter_name = images_all_txt[:,0]
all_image_paths = images_all_txt[:,1]
number_image_filters = len(all_images_filter_name)
number_images = len(all_image_paths)

image_all = np.empty(0)
image_hdu_all = np.empty(0)
image_wcs_all = np.empty(0)
for i in range(0, number_images):
	print "Opening up image: "+all_image_paths[i]
	image_all = np.append(image_all, fits.open(all_image_paths[i]))
	image_hdu_all = np.append(image_hdu_all, fits.open(all_image_paths[i])[0])
	image_wcs_all = np.append(image_wcs_all, WCS(image_hdu_all[i].header))


# Open up the photometric catalog
fitsinput = fits.open(args.input_photometry)
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
canvas=Canvas(root, height=canvasheight, width=canvaswidth)

#id_iterator = 0
#id_value = ID_numbers[id_iterator]


# Plot the EAZY SED
image = Image.open(EAZY_files+str(current_id)+"_EAZY_SED.png")

# Crop out the thumbnails
image = cropEAZY(image)

photo = resizeimage(image)
item4 = canvas.create_image(eazy_positionx, eazy_positiony, image=photo)
Label(root, text="EAZY FIT", font = "Helvetica 30").place(x=eazytext_positionx, y = eazytext_positiony)

# Plot the BEAGLE SED
new_image = Image.open(BEAGLE_files+str(current_id)+"_BEAGLE_SED.png")
new_photo = resizeimage(new_image)
item5 = canvas.create_image(beagle_positionx, beagle_positiony, image=new_photo)
Label(root, text="BEAGLE FIT", font = "Helvetica 30").place(x=beagletext_positionx, y = beagletext_positiony)
	
canvas.pack(side = TOP, expand=True, fill=BOTH)

# Plot the thumbnails
fig_photo_objects = np.empty(0, dtype = 'object')
fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, current_id, current_index, defaultstretch)


# # # # # # # # # # # #
# Flag Object Buttons

# Create the Bad Fit Flag
btn1 = Button(root, text = 'Bad Fit', bd = '5', command = badfit)
btn1.config(height = 2, width = 13, fg='red', font=('helvetica', 20))
btn1.place(x = 600, y = 870)

# Create the High Redshift Flag Button
btn1 = Button(root, text = 'High Redshift', bd = '5', command = highz)
btn1.config(height = 2, width = 13, fg='red', font=('helvetica', 20))
btn1.place(x = 800, y = 870)

# Create the Bad Data Flag
btn1 = Button(root, text = 'Bad Data', bd = '5', command = baddata)
btn1.config(height = 2, width = 13, fg='red', font=('helvetica', 20))
btn1.place(x = 1000, y = 870)


# # # # # # # # # # # #
# Move to New Object Buttons

# Create the Next Object Button
btn2 = Button(root, text = 'Next Object', bd = '5', command = nextobject)
btn2.config(height = 2, width = 20, fg='black', font=('helvetica', 20))
btn2.place(x = 1040, y = 940)

# Create the Previous Object Button
btn3 = Button(root, text = 'Previous Object', bd = '5', command = previousobject)  
btn3.config(height = 2, width = 20, fg='black', font=('helvetica', 20))
btn3.place(x = 750, y = 940)

# Create the Quit Button
btn4 = Button(root, text = 'Quit', bd = '5', command = save_destroy)  
btn4.config(height = 2, width = 20, fg='grey', font=('helvetica', 20))
btn4.place(x = 1700, y = 940)


# # # # # # # # # # # #
# Image Stretch Buttons

Label(root, text="Stretch", font = "Helvetica 20").place(x=20, y = 950)

# Create the LinearStretch Button
btn5 = Button(root, text = 'Linear', bd = '5', command = linearstretch)  
btn5.config(height = 2, width = 10, fg='grey', font=('helvetica', 20))
btn5.place(x = 100, y = 940)

# Create the SqrtStretch Button
btn6 = Button(root, text = 'Sqrt', bd = '5', command = sqrtstretch)  
btn6.config(height = 2, width = 10, fg='grey', font=('helvetica', 20))
btn6.place(x = 250, y = 940)

# Create the LinearStretch Button
btn7 = Button(root, text = 'Asinh', bd = '5', command = asinhstretch)  
btn7.config(height = 2, width = 10, fg='grey', font=('helvetica', 20))
btn7.place(x = 400, y = 940)



# Create the Object Entry Field and Button
#if (args.id_number):
#	Label(root, text="Display Object: ", font = "Helvetica 20").place(x=30, y = 940)
#	e1 = Entry(root, width = 5, font = "Helvetica 20")
#	e1.place(x = 200, y = 940)

#	btn3 = Button(root, text = 'Go', bd = '5', command = gotoobject)  
#	btn3.config(height = 2, width = 4, fg='blue', font=('helvetica', 20))
#	#btn2.pack(side = 'bottom')
#	btn3.place(x = 470, y = 940)

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
