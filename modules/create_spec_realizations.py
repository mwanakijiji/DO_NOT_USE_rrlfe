#!/usr/bin/python
'''
This module takes a list of spectra and generates normalized realizations using Gaussian Error.

This module takes a list of spectra. First it generates 100 realizations of the spectra using the error bar to move the points around
simulating Gaussian noise. Then it runs Ken Carrell's bkgrnd routine to determine the normalization and finally creates the normalized spectra.

@package create_spec_realizations
@author deleenm
@version \e \$Revision$
@date \e \$Date$

Usage: create_spec_realizations.py
'''

# -----------------------------
# Standard library dependencies
# -----------------------------
import argparse
import os
import glob
from subprocess import Popen, PIPE
## ## import test.so # experimenting with C extensions
# -------------------
# Third-party imports
# -------------------
#from astropy.io import fits ## Pylint says this is not used
from astropy.io import fits
from astropy.table import Table
#import pyfits
import sys
import numpy as np
import pandas as pd
from pathlib import *

current_dir = os.path.dirname(__file__)
print("Current dir:")
print(current_dir)

from . import *



# --------------------
# Function Definitions
# --------------------
def create_norm_spec(name_list,
                     normdir,
                     finaldir):
    '''
    Create final normalized spectra, using the output from the bkgrnd routine (which
    puts out wavelength, flux, and continuum flux, but not the actual normalized flux)

    Arguments:
        name_list: List of Realization file names (no path info)
        normdir: bkgrnd ascii files
        finaldir: The final directory for files which have completed the full normalization process.
    Returns:
       A list of final file names
    '''

    logging.info("Creating normalized spectra")

    new_name_list = list()

    for spec in name_list: # loop through spectrum realizations
        #import ipdb; ipdb.set_trace()

        # spectrum realization file name (as output by bkgrnd), with relative path info
        spec_name = os.path.join(normdir, spec)

        # astropy table containing a spectrum's 1.) wavelength, 2.) flux, 3.) background flux
        spec_tab = read_bkgrnd_spec(spec_name)
        # output file name of final, normalized spectrum, with relative path info
        new_name = os.path.join(finaldir, spec) # vestigial, from adding .smo to help Robospect pick it out
        # add to list
        new_name_list.append(new_name)
        #import ipdb; ipdb.set_trace()

        try:
            # open file to write normalized spectrum to
            outfile = open(new_name, 'w')
        except IOError:
            logging.info("File {} could not be opened!".format(new_name))
        for j in range(len(spec_tab['wavelength'])):
            # do the division to normalize and write out
            outfile.write("{} {:.4f}\n".format(spec_tab['wavelength'][j],
                                               spec_tab['flux'][j]/spec_tab['bckgrnd_flux'][j]))

        outfile.close()

    return(new_name_list)

def calc_noise(noise_level, spectrum_df):

    if (noise_level == "file"):
        # add Gaussian error to the empirical flux, taking sigma to be the
        # 'error' column in an input file; note this is an ABSOLUTE error
        noise_to_add = np.random.standard_normal(len(spectrum_df))*spectrum_df['error']
        logging.info("Injecting Gaussian noise based on error column in file.")
    elif (noise_level == "None"):
        # don't inject noise at all (note this doesn't make sense if multiple spectra are being realized)
        noise_to_add = 0
        logging.info("Injecting no noise at all")
    else:
        # noise is a set value, representing a Gaussian sigma; this is normalized,
        # so a '0.01' means 'Gaussian-distributed random number with sigma=0.01*flux_input';
        # note this 0.01 is a RELATIVE error
        noise_to_add = np.random.standard_normal(len(spectrum_df))*noise_level*spectrum_df["flux"]
        logging.info("Injecting Gaussian noise based on fixed value.")

    return noise_to_add


def generate_realizations(spec_name, outdir, spec_file_format, num, noise_level):
    '''
    Calculates a Number of Realizations of a given spectrum using Gaussian Errorbars

    Arguments:
        spec_name: The spectrum filename
        outdir: The working directory
        spec_file_format: The format of the input spectra ["fits", "ascii.no_header"]
        num: Number of realizations to generate
        noise: 'None': add no noise; 'file': take Gaussian samples of error with spread based on the error column in file
    Returns:
       A list of filenames for the realization spectra.
    '''

    logging.info("Generating spectrum realizations of " + spec_name)

    # astropy table containing an empirical spectrum's 1.) wavelength, 2.) flux, 3.) error
    # and the header of the source FITS file
    spec_tab, hdr = read_spec(spec_name, format=spec_file_format)

    #import ipdb; ipdb.set_trace()
    basename = os.path.basename(spec_name) # shave off path stem

    # generate realizations
    new_basename_list = list()

    for i in range(num):

        # make file names

        # basename of spectrum realization, ascii
        new_prefix_ascii = "{}_noise_ver_{:03d}".format(basename.split(".")[-2], i)
        suffix_ascii = basename.split(".")[-1] # could be .dat, .csv, .txt, etc.
        new_basename_ascii = new_prefix_ascii + "." + suffix_ascii
        # if we want to change to FITS intermediary files:
        new_basename_fits = "{}_noise_ver_{:03d}".format(basename.split(".fits")[0], i) + ".fits"
        # don't need path info in spec_name list; add ascii name here
        new_basename_list.append(new_basename_ascii)

        # name of spectrum realization, with path
        new_name_ascii = os.path.join(outdir, new_basename_ascii)
        new_name_fits = os.path.join(outdir, new_basename_fits)

        noise_to_add = calc_noise(noise_level=noise_level, spectrum_df=spec_tab)

        # add the noise
        new_flux = noise_to_add + spec_tab['flux']
        #print(new_flux)
        print(new_name_ascii)


        try:
            outfile = open(new_name_ascii, 'w')
        except IOError: # pragma: no cover
            logging.info("File {} could not be opened!".format(new_name_ascii))

        #import ipdb; ipdb.set_trace()

        ### write out new realization of file, in two formats:

        df_realization = Table([spec_tab["wavelength"], new_flux],
                                names=("wavelength","new_flux"))
        # set column names and write
        c1=fits.Column(name="wavelength", format='D', array=spec_tab["wavelength"])
        c2=fits.Column(name="new_flux", format='D', array=new_flux)
        t = fits.BinTableHDU.from_columns([c1, c2], header=hdr)

        ## second format: ascii, so bkgrnd can read it in
        logging.info("Writing out ascii realization file " + new_name_ascii + \
            " with noise level stdev " + str(np.std(noise_to_add)))
        for j in range(len(new_flux)):
            outfile.write("{} {:.2f}\n".format(spec_tab['wavelength'][j], new_flux[j]))
        outfile.close()

    return(new_basename_list)

def read_bkgrnd_spec(spec_name):
    '''
    Reads in ascii spectra created by bckgrnd and returns numpy arrays of wavelength, flux, bckgrnd_flux

    Arguments:
        spec_name: The spectrum filename. If Ascii file should have 3 columns: wavelength, flux, bckgrnd_flux
    Returns:
       spec_tab: A numpy Table with three columns: wavelength, flux, background flux
    '''

    logging.info("Reading ascii spectrum realization and background in " + spec_name)

    spec_tab = Table.read(spec_name, format='ascii.no_header',
                          names=['wavelength', 'flux', 'bckgrnd_flux'])

    return(spec_tab)


#####

'''
    ## note Feb. 9 2021: change to read in FITS files, but output same
    if (format == "fits"):
        # read in table info
        spec_tab = Table.read(spec_name, format='fits')
        # read in header info
        hdul = fits.open(spec_name)
        hdr = hdul[0].header
    elif (format == "ascii.no_header"):
        # this option is antiquated; we want to preserve FITS header info
        spec_tab = Table.read(spec_name, format='ascii.no_header',
                          names=['wavelength', 'flux', 'error'])
        hdr = np.nan

    return(spec_tab, hdr)
'''
#####

def read_list(input_list):
    '''
    Reads in list of spectrum names and returns a table of filenamse

    Arguments:
        input_list: a csv file with columns
            [0]: filename
            [1]: subtype (RRab, RRc)
            [2]: phase (0. to 1., -9999 for NaN)
            [3]: metallicity (if producing the calibration)
            [4]: error in metallicity (if producing the calibration)
    Returns:
       Numpy array of filenames
    '''

    logging.info("Reading in list of spectrum names to return table of filenames")

    # expects header reading
    # spectrum,subtype,phase,feh,err_feh
    input_data_arr = pd.read_csv(input_list)

    # col 0 contains the file names
    filenames_arr = input_data_arr["orig_spec_file_name"].values

    return(filenames_arr)


def read_spec(spec_name, format):
    '''
    Reads in ascii empirical spectra and returns numpy arrays of
    wavelength, flux, and error.

    Arguments:
        spec_name: The spectrum filename. If Ascii file should have
           3 columns (wavelength, flux, error no headers)
        format: "fits" or "ascii.no_header"
    Returns:
       spec_tab: A numpy Table with three columns: wavelength, flux, error
       hdr: FITS header of the input spectrum
    '''

    logging.info("Reading spectrum " + spec_name)

    if (format == "ascii.no_header"):
        # this option is antiquated; we want to preserve FITS header info
        spec_tab = Table.read(spec_name, format='ascii.no_header',
                          names=['wavelength', 'flux', 'error'])
        hdr = np.nan

    else: # pragma: no cover
        logging.info("File format unknown!!!")

    return(spec_tab, hdr)


def write_bckgrnd_input(name_list, indir, normdir):
    '''
    Create input file for the bckgrnd program

    Arguments:
        name_list: List of Realization file names (no path info)
        indir: The working directory with files to be fed into bkgrnd routine
        normdir: The output directory for normalized files
    Returns:
       A string with the background input filename; the filename itself which
       has been written out lists the input and output directories, and a
       list of the FITS files which are the spectrum realizations in the input
       directory
    '''

    logging.info("Creating input file list of spectrum realization filenames")

    #Check to see if inputfile is already there
    bckgrnd_input = os.path.join(indir, "bckgrnd_input.txt")
    if os.path.isfile(bckgrnd_input) is True:
        os.remove(bckgrnd_input)
    try:
        outfile = open(bckgrnd_input, 'w')
    except IOError:
            logging.info("File {} could not be opened!".format(bckgrnd_input))


    #Write the text file header (in_dir out_dir)
    outfile.write("{} {}\n".format(indir, normdir))
    for j in range(len(name_list)):
        outfile.write("{}\n".format(name_list[j]))
    outfile.close()

    return(bckgrnd_input)

# -------------
# Main Function
# -------------
def create_spec_realizations_main(noise_level,
                                spec_file_type,
                                num = 100,
                                  input_spec_list_dir = config_red["data_dirs"]["DIR_SRC"],
                                  input_list = config_red["data_dirs"]["DIR_SRC"] + config_red["file_names"]["LIST_SPEC_PHASE"],
                                  unnorm_empirical_spectra_dir = config_red["data_dirs"]["DIR_RAW_SPEC_DATA"],
                                  unnorm_noise_churned_spectra_dir = config_red["data_dirs"]["DIR_REZNS_SPEC"],
                                  bkgrnd_output_dir = config_red["data_dirs"]["DIR_REZNS_SPEC_NORM"],
                                  final_dir = config_red["data_dirs"]["DIR_REZNS_SPEC_NORM_FINAL"],
                                  verb=False):
    '''
    INPUTS:
    num: number of spectrum realizations to make, per empirical spectrum
    spec_file_type: file format of input spectra ["fits"/"ascii.no_header"]
    input_spec_list_dir: directory containing list of empirical spectra (## OBSOLETE? ##)
    input_list: file listing spectra we want to normalize
    unnorm_empirical_spectra_dir: directory of empirical spectra (or, if they are actually
        synthetic spectra, these are the original synthetic spectra which we will generate
        multiple realizations of)
    unnorm_noise_churned_spectra_dir: directory to contain noise-churned spectrum realizations
    bkgrnd_output_dir: directory to contain output of bkgrnd (spectra and fit continuua)
    final_dir: directory to contain normalized spectrum realizations

    OUTPUTS:
    (text files written)
    '''

    logging.info("--------------------------")
    logging.info("Making "+str(num)+" realizations of each input spectrum")

    if (num > 1) and (noise_level == "None"):
        logging.warning("Realizing multiple spectra but noise level is zero")
        input("Hit [Enter] to continue")

    # Read list of input spectra
    # input_list ALREADY SET IN DEFAULTS ## input_list = input_spec_list_dir + config_red["file_names"]["LIST_SPEC_PHASE"]
    list_arr = read_list(input_list)
    #logging.info('list_arr')
    #logging.info(list_arr)

    # Check to make sure outdir (to receive realizations of spectra) exists
    outdir = unnorm_noise_churned_spectra_dir
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    else:
        # Check to see if there are any files in the output directories (if not,
        # there is data from a previous run that will inadvertently be used later)
        # read the entries

        #import ipdb; ipdb.set_trace()

        with os.scandir(outdir) as list_of_entries1:
            counter1 = 0
            for entry1 in list_of_entries1:
                if entry1.is_file():
                    counter1 += 1
        if (counter1 != 0):
            logging.info("------------------------------")
            logging.info("Directory to write realizations not empty!!")
            logging.info(outdir)
            logging.info("------------------------------")
            input("Do what you want with those files, then hit [Enter]")

        with os.scandir(bkgrnd_output_dir) as list_of_entries2:
            counter2 = 0
            for entry2 in list_of_entries2:
                if entry2.is_file():
                    counter2 += 1
        if (counter2 != 0):
            logging.info("------------------------------")
            logging.info("Directory to write raw normalization output not empty!!")
            logging.info(bkgrnd_output_dir)
            logging.info("------------------------------")
            input("Do what you want with those files, then hit [Enter]")

        with os.scandir(final_dir) as list_of_entries3:
            counter3 = 0
            for entry3 in list_of_entries3:
                if entry3.is_file():
                    counter3 += 1
        if (counter3 != 0):
            logging.info("------------------------------")
            logging.info("Directory to write final normalization output not empty!!")
            logging.info(final_dir)
            logging.info("------------------------------")
            input("Do what you want with those files, then hit [Enter]")

    # create noise-churned realizations for each spectrum
    name_list = list() # initialize
    #import ipdb; ipdb.set_trace()
    for i in range(len(list_arr)): # make spectrum realizations and list of their filenames
        #import ipdb; ipdb.set_trace()
        name_list.extend(generate_realizations(spec_name=unnorm_empirical_spectra_dir+"/"+list_arr[i],
                                               outdir=outdir,
                                               spec_file_format=spec_file_type,
                                               num=num,
                                               noise_level=noise_level))

    # next we need to normalize the spectra; begin by creating input list of
    # spectrum realizations written in the previous step
    bkg_input_file = write_bckgrnd_input(name_list, outdir, bkgrnd_output_dir)
    logging.info("-------------------------------------------")
    logging.info('The file containing the list of spectra which will be fed ' + \
                'into the normalization routine is ' + bkg_input_file)

    # normalize each spectrum realization (smoothing parameter is set in __init__)
    bkgrnd = Popen([str(config_red["data_dirs"]["DIR_BIN"]) + "bkgrnd", "--smooth "+str(config_red["reduc_params"]["SMOOTH"]),
                    "--sismoo 1", "--no-plot", "{}".format(bkg_input_file)], stdout=PIPE, stderr=PIPE)
    (out, err) = bkgrnd.communicate() # returns tuple (stdout, stderr)

    if verb == True: ## ## decode messages (are they used later? why take this step?)
        logging.info(out.decode("utf-8"))
        logging.info(err.decode("utf-8"))

    # read in input files, normalize them, write out, and return list of those filenames
    final_list = create_norm_spec(name_list, bkgrnd_output_dir, final_dir)

    logging.info("-------------------------------------------")
    logging.info("Wrote realizations of original spectra to directory")
    logging.info(outdir)
    logging.info("-------------------------------------------")
    logging.info("Wrote raw normalization output to directory")
    logging.info(bkgrnd_output_dir)
    logging.info("-------------------------------------------")
    logging.info("Wrote final normalized spectra to directory")
    logging.info(final_dir)
