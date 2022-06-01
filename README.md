# JADESView

```
% python JADESView.py -h
usage: JADESView.py [-h] [-input INPUT] [-id ID_NUMBER]
                    [-idlist ID_NUMBER_LIST]

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT          JADESView Input File?
  -id ID_NUMBER         ID Number?
  -idlist ID_NUMBER_LIST
                        List of ID Numbers?
  -idarglist IDARGLIST  Command line argument list of objects
```


This tool was designed to allow users to look at EAZY and BEAGLE fits alongside
thumbnails for for objects from JADES photometric catalogs. The user can either download
the photometric catalogs and the EAZY and BEAGLE plots from this website:

[https://fenrir.as.arizona.edu/jades/data/](https://fenrir.as.arizona.edu/jades/data/)

(The log-in details and the catalog descriptions are found [on confluence here](https://issues.cosmos.esa.int/jwst-nirspecwiki/display/WGs/Step+04+-+Photometric+redshifts+and+derived+information).)

Or they can point to the catalogs and files directly on fenrir. The paths to the EAZY and BEAGLE 
files on fenrir, as well as the username and password, are provided on the confluence page (NOTE: 
you must use https). It is recommended that you have a pretty speedy internet connection, since this 
will add some time (in my tests, around 1 second per object) for fetching the images from the 
server instead of on your local machine, but it saves having to download many man GBs of png files.

The tool requires numpy, matplotlib, tkinter, and astropy installations, and is written
using Python 3.0 (but will work under Python 2 as well). [You can learn more about the
packages that are required here](https://github.com/kevinhainline/JADESView#installation). 

The user will have to modify the file `JADESView_input_file.dat` to point to these files
and also to specify some default parameters:

```
input_photometry       /Path/to/Photometric_Catalog.fits
image_list             /Path/to/image_list.dat 
EAZY_files             /Path/to/EAZY_output_plots/
EAZY_results           /Path/to/EAZY_results_summary.fits
BEAGLE_files           /Path/to/BEAGLE_output_plots/
BEAGLE_results         /Path/to/Beagle_VAC.fits
output_flags_file      Object_Flags.fits
output_notes_file      Object_Notes.txt
canvasheight           1000
canvaswidth            2000
defaultstretch         AsinhStretch
ra_dec_size_value      2.0
fenrir_username        fenrir_username
fenrir_password        fenrir_password
```

In this file, do not modify the first column, but replace the values in the second column
with the desired input. The file specified by `input_photometry` is the desired photometric 
catalog created from the mosaics. The `image_list` is a file with two columns:

```
NRC_F090W /Path/to/Mosaic/File/F090W.fits
NRC_F115W /Path/to/Mosaic/File/F115W.fits
NRC_F150W /Path/to/Mosaic/File/F150W.fits
NRC_F200W /Path/to/Mosaic/File/F200W.fits
NRC_F277W /Path/to/Mosaic/File/F277W.fits
NRC_F335M /Path/to/Mosaic/File/F335M.fits
NRC_F356W /Path/to/Mosaic/File/F356W.fits
NRC_F410M /Path/to/Mosaic/File/F410M.fits
NRC_F444W /Path/to/Mosaic/File/F444W.fits
```
etc etc. This allows the program to find the individual mosaic images for the thumbnails.
Currently, it is not feasible to continually fetch the mosaic images from a given server, so
these will need to be downloaded. 

The `EAZY_files` and `BEAGLE_files` point to the folders on your local machine, or on fenrir, 
that contain the EAZY and BEAGLE output plots from the above linked repository. Similarly, if the
user specifies an `EAZY_results` (the `EAZY_results_summary.fits` file that is also found on
fenrir) and `BEAGLE_results` (the `*_BEAGLE_VAC.fits` file that is also found on fenrir) files, 
which contain the output photometric redshifts, these values will be printed to the window for
each object. 

Note that the current version will allow the user to change the canvaswidth (default: 2000 pixels), 
and hopefully everything will  scale so that you can use the code on smaller monitors. 

The program will produce output files that may be useful. The first is a fits file with 
objects flagged (currently there are three flags, bad data, bad fit, and high-redshift object,
but this can be expanded, if anyone has any suggestions), and the user specifies the name
of this output file with the `output_flags_file` name. The user can also write notes 
within the program on individual objects which will be put into a text file named based 
on the entry in the `Object_Notes.txt` file. Currently, the program was written to be
1000x2000 pixels (and all of the button and plot placement supports that), so I wouldn't
change the canvasheight or canvaswidth parameters yet. The user can change the default
stretch on the thumbnails with the `defaultstretch` entry (current options are: `AsinhStretch`, 
`LinearStretch`, and `LogStretch`). There is a button that the user can press which will add a 
crosshair on each of the thumbnails at the location of the current target. Finally, the default size 
in arcseconds of the thumbnails is given by `ra_dec_size_value.` Both the stretch and the thumbnail 
size can be changed within the program on the fly. 

The program is run by specifying an ID, an ID list (as a text file), or an ID list on the
command line:

```
python JADESView.py -id 1234
python JADESView.py -idlist list-of-IDs.dat
python JADESView.py -idarglist '5 6 7 8 9 10' 
```

If no IDs are provided, the code will start at the first object in the `input_photometry` ID list.  

Also, the user can specify a different input file using the `-input` flag, otherwise
it will default to looking for the file `JADESView_input_file.dat.`

Within the program, the user can look at the EAZY and BEAGLE SED fits, chi-square surfaces,
and posteriors, alongside the thumbnails (where the flux is not given as -9999), with
the SNR in each filter provided at the bottom of each panel. Below the thumbnails, the
user can change the size of the thumbnail, as well as the stretch. In the center, along
the bottom, are buttons to flag individual objects, as well as advance to the next or return
to the previous object (if using a list of targets, this will move along the list, otherwise
it will advance numerically to the next or previous object). If the user has not specified
a list of targets, there will also be an option to go to any object with a given ID. 
Finally, there's a field for writing notes, which will be saved until the program is 
quit, after which the notes will be written to a file. In addition, the user can save the
contents of the canvas (the EAZY, BEAGLE fits, and thumbnails only, not the buttons)
to an image file ('XXXX_JADESView.png' where XXXX is the Object ID) with the Save Canvas 
button. 

Note: The program will overwrite the output files if they are not renamed between runs. 

## Installation

While JADESView works under python3, you will need to use matplotlib v2.2.5. I recommend creating
a conda environment using the attached `environment.yml` file:

```
% conda env create -f environment.yml
```

This will create a new conda environment, `jadesview`, from which you can run `JADESView` 
and `JADESView_webget`.

...or you could use the `requirements.txt` file to install these packages using `pip`. 