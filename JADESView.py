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

EAZY_files = '/Volumes/KNH_EXTERNAL/December_2019_DataChallenge/EAZY_run_1_13_20/F200W_plots_kron_f80/'
BEAGLE_files = '/Volumes/KNH_EXTERNAL/December_2019_DataChallenge/BEAGLE_output_1_13_20/F200W_kron_f80/'

def resizeimage(image):
	basewidth = 1000
	wpercent = (basewidth / float(image.size[0]))
	hsize = int((float(image.size[1]) * float(wpercent)))
	image = image.resize((basewidth, hsize), PIL.Image.ANTIALIAS)
	photo = ImageTk.PhotoImage(image)
	return photo
	

def highz():
	print "This is a High-Redshift Candidate"


def nextobject():
	#print "Moving to the next object."
	global id_iterator
	global e2
	global ID_numbers
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

	notes_values[id_iterator] = e2.get()
	e2.delete(0,END)

	canvas.delete(item4)
	canvas.delete(item5)
	if (id_iterator < len(ID_numbers)-1):
		id_iterator = id_iterator+1
	else:
		id_iterator = len(ID_numbers)-1
	
	id_value = ID_numbers[id_iterator]
	e2.insert(0, notes_values[id_iterator])

	#print id
	image = Image.open(EAZY_files+str(id_value)+"_EAZY_SED.png")
	image = cropEAZY(image)
	photo = resizeimage(image)
	item4 = canvas.create_image(eazy_positionx, eazy_positiony, image=photo)
	
	new_image = Image.open(BEAGLE_files+str(id_value)+"_BEAGLE_SED.png")
	new_photo = resizeimage(new_image)
	item5 = canvas.create_image(beagle_positionx, beagle_positiony, image=new_photo)

	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, id_value, defaultstretch)


def previousobject():
	#print "Moving to the next object."
	
	global id_iterator
	global e2
	global ID_numbers
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

	notes_values[id_iterator] = e2.get()
	e2.delete(0,END)
	
	canvas.delete(item4)
	canvas.delete(item5)
	if (id_iterator >= 0):
		id_iterator = id_iterator-1
	else:
		id_iterator = 0
	
	id_value = ID_numbers[id_iterator]
	e2.insert(0, notes_values[id_iterator])
	
	#print id
	image = Image.open(EAZY_files+str(id_value)+"_EAZY_SED.png")
	image = cropEAZY(image)
	photo = resizeimage(image)
	item4 = canvas.create_image(eazy_positionx, eazy_positiony, image=photo)
	
	new_image = Image.open(BEAGLE_files+str(id_value)+"_BEAGLE_SED.png")
	new_photo = resizeimage(new_image)
	item5 = canvas.create_image(beagle_positionx, beagle_positiony, image=new_photo)

	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, id_value, defaultstretch)


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
	global id_iterator
	global ID_numbers
	global canvas   
	global fig_photo_objects
	global defaultstretch

	defaultstretch = 'LinearStretch'
	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, ID_numbers[id_iterator], defaultstretch)

def sqrtstretch():
	global id_iterator
	global ID_numbers
	global canvas   
	global fig_photo_objects
	global defaultstretch

	defaultstretch = 'SqrtStretch'
	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, ID_numbers[id_iterator], defaultstretch)

def asinhstretch():
	global id_iterator
	global ID_numbers
	global canvas   
	global fig_photo_objects
	global defaultstretch

	defaultstretch = 'AsinhStretch'
	fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, ID_numbers[id_iterator], defaultstretch)


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

def create_thumbnails(canvas, fig_photo_objects, id_value, stretch):
	global ID_values
	global RA_values
	global DEC_values
	global decsize_value
	global rasize_value
	global image_all
	global image_hdu_all
	global image_wcs_all
	global image_flux_value_err_cat
	global all_images_filter_name
	global number_images
	global SNR_values
	
	# Let's associate the selected object with it's RA and DEC		
	# Create the object thumbnails. 
	idx_cat = np.where(ID_values == id_value)[0]
	objRA = RA_values[idx_cat][0]
	objDEC = DEC_values[idx_cat][0]
			
	cosdec_center = math.cos(objDEC * 3.141593 / 180.0)
		
	# Set the position of the object
	position = SkyCoord(str(objRA)+'d '+str(objDEC)+'d', frame='fk5')
	size = u.Quantity((decsize_value, rasize_value), u.arcsec)
	
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
			ax3.text(0.5, 0.95, all_images_filter_name[i].split('_')[1], transform=ax3.transAxes, fontsize=10, fontweight='bold', ha='center', va='top', color = 'white')
			ax3.text(0.95, 0.05, 'SNR = '+str(round(SNR_values[idx_cat, i],2)), transform=ax3.transAxes, fontsize=7, fontweight='bold', horizontalalignment='right', color = 'white')
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
							
			# Keep this handle alive, or else figure will disappear
			if (i <= 5):
				fig_x, fig_y = 20+(175*i), 500
				#Label(root, text=all_images_filter_name[i].split('_')[1], font = "Helvetica 20").place(x=fig_x, y = fig_y)
			if ((i > 5) & (i <= 11)):
				fig_x, fig_y = 20+(175*(i-6)), 675
				#Label(root, text=all_images_filter_name[i].split('_')[1], font = "Helvetica 20").place(x=fig_x, y = fig_y)
			if ((i > 11) & (i <= 17)):
				fig_x, fig_y = 20+(175*(i-12)), 900
				#Label(root, text=all_images_filter_name[i].split('_')[1], font = "Helvetica 20").place(x=fig_x, y = fig_y)
				
			fig_photo_objects = np.append(fig_photo_objects, draw_figure(canvas, fig, loc=(fig_x, fig_y)))
			plt.close('all')

	return fig_photo_objects



def save_destroy():
	global id_iterator
	
	notes_values[id_iterator] = e2.get()
	
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

# Open up the filters files
all_filters_file_name = '/Users/knh/Desktop/NIRCam/photometric_redshifts/PhotoZReady/JADES_All_Filters.dat'
filters_all = np.loadtxt(all_filters_file_name, dtype='str')
all_filters = filters_all[:,0]
all_filter_numbers = filters_all[:,1]
all_header_filters = filters_all[:,2]
all_filter_waves = filters_all[:,3]
all_filter_waves = all_filter_waves.astype(np.float)
all_filter_bw = filters_all[:,4]
all_filter_bw = all_filter_bw.astype(np.float)
total_all_jades_filters = len(all_filters)

# Now, what are the filters that we're going to actually use?
filter_file_name = args.filters
filters_file = np.loadtxt(filter_file_name, dtype='str')
filters = filters_file
number_filters = len(filters)

filter_numbers = np.zeros(number_filters, dtype = int)
header_filters = np.zeros(number_filters, dtype = '|S30')
filter_waves = np.zeros(number_filters)
filter_bw = np.zeros(number_filters)
for n in range(0, number_filters):
	filter_index = np.where(all_filters == filters[n])[0][0]
	header_filters[n] = all_header_filters[filter_index]
	filter_numbers[n] = all_filter_numbers[filter_index]
	filter_waves[n] = all_filter_waves[filter_index]
	filter_bw[n] = all_filter_bw[filter_index]

all_images_file_name = args.image_list
images_all_txt = np.loadtxt(all_images_file_name, dtype='str')
all_images_filter_name = images_all_txt[:,0]
all_image_paths = images_all_txt[:,1]
number_image_filters = len(all_images_filter_name)
number_images = len(all_image_paths)


# Open up the photometric catalog
fitsinput = fits.open(args.input_photometry)
ID_values = fitsinput[1].data['ID']
RA_values = fitsinput[1].data['RA']
DEC_values = fitsinput[1].data['DEC']
number_objects = len(ID_values)

image_flux_value_cat = np.zeros([number_objects, number_image_filters])
image_flux_value_err_cat = np.zeros([number_objects, number_image_filters])
SNR_values = np.zeros([number_objects, number_filters])

for j in range(0, number_image_filters):
	image_flux_value_cat[:,j] = fitsinput[1].data[all_images_filter_name[j]]
	image_flux_value_err_cat[:,j] = fitsinput[1].data[all_images_filter_name[j]+'_err']
	SNR_values[:,j] = image_flux_value_cat[:,j] / image_flux_value_err_cat[:,j]


rasize_value = 2.0
decsize_value = 2.0

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
	#image_all = np.append(image_all, fits.open(all_image_paths[i])[0].data)
	image_hdu_all = np.append(image_hdu_all, fits.open(all_image_paths[i])[0])
	image_wcs_all = np.append(image_wcs_all, WCS(image_hdu_all[i].header))


if (args.id_number):
	ID_numbers = np.zeros(1).astype(int)
	ID_numbers[0] = args.id_number
	number_input_objects = len(ID_numbers)
	if (args.id_number_list):
		print "You can't specify an individual ID and a list, ignoring the list."

if not (args.id_number):
	if not (args.id_number_list):
		sys.exit("Need to specify an ID number or a list of ID numbers!")

if (args.id_number_list):
	ID_input_file = np.loadtxt(args.id_number_list)
	if (len(ID_input_file.shape) > 1):
		ID_numbers = ID_input_file[:,0].astype(int)
	else:
		ID_numbers = ID_input_file.astype(int)
	#ID_numbers = ID_input_file
	number_input_objects = len(ID_numbers)

notes_values = np.array([''], dtype = 'object')
for x in range(0, number_input_objects-1):
	notes_values = np.append(notes_values, ['']) 

# Start by creating the GUI as root
root=Tk()
root.wm_title("JADESView")

id_iterator = 0
id_value = ID_numbers[id_iterator]

# Plot the EAZY SED
image = Image.open(EAZY_files+str(id_value)+"_EAZY_SED.png")

eazy_positionx, eazy_positiony = 500, 230
eazytext_positionx, eazytext_positiony = 350, 70
beagle_positionx, beagle_positiony = 1500, 350
beagletext_positionx, beagletext_positiony = 1300, 70

# For the future, for removing 
image = cropEAZY(image)
canvas=Canvas(root, height=1000, width=2000)
photo = resizeimage(image)
item4 = canvas.create_image(eazy_positionx, eazy_positiony, image=photo)
Label(root, text="EAZY FIT", font = "Helvetica 30").place(x=eazytext_positionx, y = eazytext_positiony)

# Plot the BEAGLE SED
new_image = Image.open(BEAGLE_files+str(id_value)+"_BEAGLE_SED.png")
new_photo = resizeimage(new_image)
item5 = canvas.create_image(beagle_positionx, beagle_positiony, image=new_photo)
Label(root, text="BEAGLE FIT", font = "Helvetica 30").place(x=beagletext_positionx, y = beagletext_positiony)
	
canvas.pack(side = TOP, expand=True, fill=BOTH)

# Plot the thumbnails
defaultstretch = 'AsinhStretch'
fig_photo_objects = np.empty(0, dtype = 'object')
fig_photo_objects = create_thumbnails(canvas, fig_photo_objects, id_value, defaultstretch)



# Create the High Redshift Flag Button
btn1 = Button(root, text = 'Flag As High Redshift', bd = '5', command = highz)
btn1.config(height = 2, width = 20, fg='red', font=('helvetica', 20))
btn1.place(x = 900, y = 870)

# Create the Next Object Button
if (args.id_number):
	btn2 = Button(root, text = 'Next Object', bd = '5', command = nextobject)
	btn2.config(height = 2, width = 20, fg='grey', font=('helvetica', 20), state=DISABLED)
	btn2.place(x = 1040, y = 940)
else:
	btn2 = Button(root, text = 'Next Object', bd = '5', command = nextobject)
	btn2.config(height = 2, width = 20, fg='black', font=('helvetica', 20))
	btn2.place(x = 1040, y = 940)

# Create the Previous Object Button
if (args.id_number):
	btn3 = Button(root, text = 'Previous Object', bd = '5', command = previousobject)  
	btn3.config(height = 2, width = 20, fg='black', font=('helvetica', 20), state=DISABLED)
	btn3.place(x = 750, y = 940)
else:
	btn3 = Button(root, text = 'Previous Object', bd = '5', command = previousobject)  
	btn3.config(height = 2, width = 20, fg='black', font=('helvetica', 20))
	btn3.place(x = 750, y = 940)

# Create the Quit Button
btn3 = Button(root, text = 'Quit', bd = '5', command = save_destroy)  
btn3.config(height = 2, width = 20, fg='grey', font=('helvetica', 20))
btn3.place(x = 1700, y = 940)

Label(root, text="Stretch", font = "Helvetica 20").place(x=20, y = 950)

# Create the LinearStretch Button
btn3 = Button(root, text = 'Linear', bd = '5', command = linearstretch)  
btn3.config(height = 2, width = 10, fg='grey', font=('helvetica', 20))
btn3.place(x = 100, y = 940)

# Create the SqrtStretch Button
btn3 = Button(root, text = 'Sqrt', bd = '5', command = sqrtstretch)  
btn3.config(height = 2, width = 10, fg='grey', font=('helvetica', 20))
btn3.place(x = 250, y = 940)

# Create the LinearStretch Button
btn3 = Button(root, text = 'Asinh', bd = '5', command = asinhstretch)  
btn3.config(height = 2, width = 10, fg='grey', font=('helvetica', 20))
btn3.place(x = 400, y = 940)



# Create the Object Entry Field and Button
if (args.id_number):
	Label(root, text="Display Object: ", font = "Helvetica 20").place(x=30, y = 940)
	e1 = Entry(root, width = 5, font = "Helvetica 20")
	e1.place(x = 200, y = 940)

	btn3 = Button(root, text = 'Go', bd = '5', command = gotoobject)  
	btn3.config(height = 2, width = 4, fg='blue', font=('helvetica', 20))
	#btn2.pack(side = 'bottom')
	btn3.place(x = 470, y = 940)

# Create the Notes Field
Label(root, text="Notes", font = "Helvetica 20").place(x=1200, y = 870)
e2 = Entry(root, width = 50, font = "Helvetica 20")
e2.place(x = 1300, y = 870)
e2.insert(0, notes_values[id_iterator])




root.mainloop()
