# JADESView

```
% python JADESView.py -h
usage: JADESView.py [-h] -jv_params JADESVIEW_INPUT_FILE [-zlist Z_NUMBER_LIST] 
	[-id ID_NUMBER] [-idlist ID_NUMBER_LIST] [-idarglist IDARGLIST]

options:
  -h, --help            show this help message and exit
  -jv_params JADESVIEW_INPUT_FILE, --jv_params JADESVIEW_INPUT_FILE
                        JADESView input parameter file
  -zlist Z_NUMBER_LIST, --z_number_list Z_NUMBER_LIST
                        List of ID numbers and redshfits for a subsample
  -id ID_NUMBER, --id_number ID_NUMBER
                        ID Number
  -idlist ID_NUMBER_LIST, --id_number_list ID_NUMBER_LIST
                        List of ID Numbers
  -idarglist IDARGLIST  Command line argument list of objects
```


This tool was designed to allow users to look at EAZY fits alongside
thumbnails for for objects from JADES photometric catalogs. The user should download
the latest photometric catalog from here:

[https://archive.stsci.edu/hlsp/jades](https://archive.stsci.edu/hlsp/jades)

JADESView will require the installation of [EAZY-py](https://github.com/gbrammer/eazy-py), 
as it serves as a wrapper and visualizer for running this photometric redshift code on
sources within the JADES catalog. In addition, the tool requires a number of other
python packages, and is written using Python 3.12+. It is recommended that you use
either [conda](https://www.anaconda.com/) or 
[micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) 
to set up an environment for running JADESView. [You can learn more about the
packages that are required here](https://github.com/kevinhainline/JADESView#installation). 

The user will have to modify the file `JADESView_input_file.dat` to point to these files
and also to specify some default parameters:

```
input_photometry         /Path/to/photometric/catalog
convolved                False
aperture                 CIRC1
uncertainty              _e
filters                  JADES_filters_no_WFC3_v1.0.0.dat
min_rel_err              0.05
eazy_ncores              8
ancillary_files          JADESView_files
zeropoint                zphot.zeropoint
templates                JADES_fsps.param
tempfilt                 JADES_GS_v100_tempfilt.pkl 
image_list               JADES_v100_mosaic_files.dat
plot_thumbnails          True
ra_dec_size_value        2.0
asada_cgm                False
fix_zspec                False
fitsmap_link             https://jades.idies.jhu.edu/?ra=RA_VALUE&dec=DEC_VALUE&zoom=10
```

In this file, do not modify the first column, but replace the values in the second column
with the desired input. The file specified by `input_photometry` is the JADES hlsp
catalog. 

`JADESView_input_file.dat` file includes the option to change the photometry being fit: convolved, and
which apertures, and uncertainties within the catalog. The `filters` file is the list
of filters to be fit. For instance, here is a large sample of filters that could be
fit within JADES:

```
HST_F435W
HST_F606W
NRC_F070W
HST_F775W
HST_F814W
HST_F850LP
NRC_F090W
NRC_F115W
NRC_F150W
NRC_F162M
NRC_F182M
NRC_F200W
NRC_F210M
NRC_F250M
NRC_F277W
NRC_F300M
NRC_F335M
NRC_F356W
NRC_F410M
NRC_F430M
NRC_F444W
NRC_F460M
NRC_F480M
```

If there is no flux for this filter for a given objet in the catalog, this will be set
to -9999 +/- -9999 and ignored by EAZY. Note that the filters being fit with EAZY do not
necessarily have to match the image filters - you can include as many images of the 
region, providing you're confident that their astrometic solution matches the JADES
survey. You can also point to a `zeropoint` file, which is generaly of the form `zphot.zeropoint`
and lists the photometric offsets that will be applied when running EAZY. This doesn't
have to be used, and can be set to a nonexistent file, where it will be ignored.

One thing that is important for JADESView is that, to speed up the fit, it 
will look for a `tempfilt` file, generally in pickle format. This file is generated for a
specific filter set, template set, EAZY Z_MIN, Z_MAX, Z_STEP, and CGM prescription. If
you change any of these parameters, and try to use a tempfilt file that does not match,
the code will exit with an error message. The length of time it takes to create the
tempfilt array depends on the number of templates, filters, and redshift grid parameters.

You can also plot thumbnails underneath the SED and chi-square surface, by
setting `plot_thumbnails` to `True`. If you do, you should point to the various 
mosaic values with the `image_list`. The code looks for three columns:

```
HST_F435W 0 /Path/to/Mosaic/File/F435W.fits
NRC_F090W 1 /Path/to/Mosaic/File/F090W.fits
NRC_F115W 1 /Path/to/Mosaic/File/F115W.fits
NRC_F150W 1 /Path/to/Mosaic/File/F150W.fits
NRC_F200W 1 /Path/to/Mosaic/File/F200W.fits
NRC_F277W 1 /Path/to/Mosaic/File/F277W.fits
NRC_F335M 1 /Path/to/Mosaic/File/F335M.fits
NRC_F356W 1 /Path/to/Mosaic/File/F356W.fits
NRC_F410M 1 /Path/to/Mosaic/File/F410M.fits
NRC_F444W 1 /Path/to/Mosaic/File/F444W.fits
```
The first column is the name of the filter (the names that are supported are in column 1
of `JADESView_files/JADES_All_Filters_EAZYpy.dat`. The second column is the hdu of the
data, and the third column is the path. The default size of the thumbnails, in arcseconds,
is given by `ra_dec_size_value.` The thumbnail size can be changed within the program on 
the fly. 

The program is run by specifying the `JADESView_input_file.dat` and an ID, an ID list
 (as a text file), or an ID list on the command line:

```
python JADESView.py -jv_params JADESView_input_file.dat -id 1234
python JADESView.py -jv_params JADESView_input_file.dat -idlist list-of-IDs.dat
python JADESView.py -jv_params JADESView_input_file.dat -idarglist '5 6 7 8 9 10' 
```

Within the program, the user can look at the EAZY SED fit, chi-square surface,
alongside the thumbnails (if specified), with the SNR in each filter provided at the 
bottom of each panel, if the flux can be gleaned from the photometric catalog. 

At the bottom of the screen is a series of buttons and text entries. On the top row 
are object ID (you can only insert the ones in your input list), X and Y min and max
values for the SED plot, a text box where you can enter another alternate redshift to
look at the best EAZY fit at that redshift, and then a previous and next button for 
cycling through your input list (you can also use PgUp and PgDown on your keyboard
to do this). On the right you can set the thumbnail size (in arcseconds),
you can click to change the chi-square surface to a P(z) plot, and you can set the
maximum Y value for the chisq surface.

On the bottom row, the leftmost button updates the plots based on your new entries (although,
at any point, you can also press "return" or "enter" on your keyboard to do the same). Next
to that is a button that resets the axes of the plots.

On the right, in the bottom row, are a button that points to the `fitsmap` location of
the object, given a fitsmap link provided in `JADESView_input_file.dat`. Notice in the
example:

```
fitsmap_link             https://jades.idies.jhu.edu/?ra=RA_VALUE&dec=DEC_VALUE&zoom=10
```

The words `RA_VALUE` and `DEC_VALUE`. These are what the code is looking for, and so
if the fitsmap button isn't working, check your file. 

The next button, `Save Thumbnails`, allow the user to create an output set of thumbnails as fits file 
cutouts in a folder labelled after the ID for the object. Similarly, `Save EAZY SED + Chisq`, 
will, at the user specified redshift or the best-fit redshift depending on what's selected,
will produce a pair of text files that have the EAZY SED (wavelength in microns and flux
in nJy) and the chisq surface (redshift and chisq). The third button `Save as PNG`, will
save the whole GUI to a PNG, also labelled with the ID. Try it out! 

Finally, as you might expect, `Quit` quits the program. 

## Example

If you want to see the power of the code, you can download the "NIRCam Photometry (GOODS-S-Deep v2.0)"
file from [the public JADES page](https://archive.stsci.edu/hlsp/jades). This file is named
`hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits`, and you can point to
this file and its location within `JADESView_input_file.dat`. Now, you can use the files
`JADESView_files/JADES_filters_no_WFC3_v1.0.0.dat` for the filter selection, and 
`JADESView_files/templates/JADES_fsps.dat`, which correspond to `JADES_GS_v100_tempfilt.pkl`,
so your `JADESView_input_file.dat` would look like this:

```
input_photometry         hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits
convolved                False
aperture                 CIRC1
uncertainty              _e
filters                  JADES_filters_no_WFC3_v1.0.0.dat
min_rel_err              0.05
eazy_ncores              8
ancillary_files          JADESView_files
zeropoint                zphot.zeropoint
templates                JADES_fsps.param
tempfilt                 JADES_GS_v100_tempfilt.pkl 
image_list               JADES_v100_mosaic_files.dat
plot_thumbnails          False
ra_dec_size_value        2.0
asada_cgm                False
fix_zspec                False
fitsmap_link             https://jades.idies.jhu.edu/?ra=RA_VALUE&dec=DEC_VALUE&zoom=9
```

And now you can run:

```
% python -W ignore JADESView.py -jv_params ../JADESView_test/JADESView_input_file.dat -idarglist '128771, 130158, 183348'
```

To look at the SEDs. If you download the very large mosaic files on STScI JADES MAST site,
you can point to those as indicated in `JADESView_files/JADES_v100_mosaic_files.dat` and
those will be plotted beneath. 

## Installation

While JADESView works under python3, you will need to use matplotlib v2.2.5. I recommend creating
a conda environment using the attached `environment.yml` file:

```
% conda env create -f environment.yml
```

This will create a new conda environment, `jadesview`, from which you can run `JADESView`.

You'll also want to install eazy-py in this environment:
```
% conda activate jadesview
% pip install eazy
% pip install git+https://github.com/karllark/dust_attenuation.git
```

## Older Version

The original version of the code can be found in the legacy-python2 branch, but it is not
recommended to use this version! 


