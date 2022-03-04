import matplotlib
matplotlib.use('Agg')

import sys, os
from configparser import ConfigParser, ExtendedInterpolation
import pandas as pd
import astropy
import pytest

current_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(current_dir, "../"))
sys.path.insert(0, target_dir)

# import more things with changed system path
from modules import *
from modules import create_spec_realizations
from conf import *
import numpy as np

# configuration data for reduction
config_red = ConfigParser(interpolation=ExtendedInterpolation()) # for parsing values in .init file
# config for reduction to find a, b, c, d
config_red.read(os.path.join(os.path.dirname(__file__), '../conf', 'config_red.ini'))

'''
# check if the directory-making function works
def test_make_dirs():

    # call function to make directories
    make_dirs() # this leads to permissions errors in online build

    # do all the directories exist now?
    for vals in config["data_dirs"]:
        abs_path_name = str(config["data_dirs"][vals])
        assert os.path.exists(abs_path_name)
'''

#def test_create_spec_realizations_main():
    ## ## CONTINUE HERE



def test_calc_noise():

    array_fake_file = np.random.random(10000)
    df_fake_file = pd.DataFrame(array_fake_file, columns=["error"])
    df_fake_file["flux"] = 2.

    noise_to_add_file = create_spec_realizations.calc_noise(noise_level = "file",
                                                            spectrum_df=df_fake_file)

    noise_to_add_none = create_spec_realizations.calc_noise(noise_level = "None",
                                                            spectrum_df=df_fake_file)

    noise_level_choice = 0.05
    noise_to_add_gauss = create_spec_realizations.calc_noise(noise_level = noise_level_choice,
                                                            spectrum_df=df_fake_file)


    # if we divide the error spectrum column back out, is the noise normally
    # distributed with stdev of 1.?
    spec_flux_file_divded = np.divide(noise_to_add_file,df_fake_file["error"])
    assert np.round(np.std(spec_flux_file_divded),1) == 1.00

    # is no noise really zero
    assert noise_to_add_none == 0

    # is Gaussian noise normalized if we divide out flux level and normalization constant?
    spec_flux_gauss_divded = np.divide(noise_to_add_gauss,noise_level_choice*df_fake_file["flux"])
    assert np.round(np.std(spec_flux_gauss_divded),1) == 1.00


def test_read_spec():

    # ascii format
    spec_name_ascii = config_red["data_dirs"]["TEST_DIR_SRC"] + "raw_spec/575020m05.smo"
    test_spec_tab_ascii, test_hdr_ascii = create_spec_realizations.read_spec(spec_name=spec_name_ascii, format="ascii.no_header")

    # for ascii data, there should be 3 columns of floats, and NO header
    assert np.isfinite(test_hdr_ascii) == False
    assert len(test_spec_tab_ascii.colnames) == 3 # 3 columns
    assert isinstance(test_spec_tab_ascii["wavelength"][0],np.float64)
    assert isinstance(test_spec_tab_ascii["flux"][0],np.float64)
    assert isinstance(test_spec_tab_ascii["error"][0],np.float64)


def test_generate_realizations():

    # use some test spectra
    abs_stem_src = config_red["data_dirs"]["TEST_DIR_SRC"] + "raw_spec/"
    abs_stem_bin = config_red["data_dirs"]["TEST_DIR_BIN"]

    # set fractional noise level for these tests
    noise_choice = 0.01

    # test on ascii files
    test_spec_list_ascii = [
                            abs_stem_src+"575020m05.smo",
                            abs_stem_src+"575020m10.smo",
                            abs_stem_src+"575020m15.smo"
                            ]
    # expected names
    # (note these should just be basenames)
    expected_filenames_ascii = [
                                "575020m05_noise_ver_000.smo",
                                "575020m05_noise_ver_001.smo",
                                "575020m10_noise_ver_000.smo",
                                "575020m10_noise_ver_001.smo",
                                "575020m15_noise_ver_000.smo",
                                "575020m15_noise_ver_001.smo"
                                ]
    returned_filenames_ascii = []
    for spec_num in range(0,len(test_spec_list_ascii)):
        return_names_one_spec = create_spec_realizations.generate_realizations(spec_name=test_spec_list_ascii[spec_num],
                                               outdir=abs_stem_bin,
                                               spec_file_format="ascii.no_header",
                                               num=2,
                                               noise_level=noise_choice)
        returned_filenames_ascii.extend(return_names_one_spec)

    # flatten lists
    #returned_filenames_fits = [item for sublist in t for item in sublist]
    #returned_filenames_ascii = [item for sublist in t for item in sublist]

    # check if elements in list 1 are in list 2
    #result_fits =  all(elem in returned_filenames_fits for elem in expected_filenames_fits)
    result_ascii =  all(elem in returned_filenames_ascii for elem in expected_filenames_ascii)

    # test: are the file names right?
    #assert result_fits
    assert result_ascii

    # test: are the original spectra divided by the realizations equivalent
    # to 1 + noise residuals? (this is just on ascii for now)

    # read in original spectra

    # sorting is critical here, to keep spectra and their realizations organized
    test_spec_list_ascii.sort() # in-place
    expected_filenames_ascii.sort()

    orig_spec_0 = pd.read_csv(test_spec_list_ascii[0], names=["wavel", "abs_flux", "error"], delim_whitespace=True)
    #orig_spec_1 = pd.read_csv(test_spec_list_ascii[1], names=["wavel", "abs_flux", "error"], delim_whitespace=True)

    # read in the realizations based off of the original spectra
    # realizations of orig spec 0
    realzn_spec_0_0 = pd.read_csv(abs_stem_bin + expected_filenames_ascii[0], names=["wavel", "abs_flux"], delim_whitespace=True)
    realzn_spec_0_1 = pd.read_csv(abs_stem_bin + expected_filenames_ascii[1], names=["wavel", "abs_flux"], delim_whitespace=True)
    # realizations of orig spec 1
    #realzn_spec_1_0 = pd.read_csv(abs_stem_bin + expected_filenames_ascii[2], names=["wavel", "abs_flux"], delim_whitespace=True)
    #realzn_spec_1_1 = pd.read_csv(abs_stem_bin + expected_filenames_ascii[3], names=["wavel", "abs_flux"], delim_whitespace=True)

    # do the division
    div_spec_0_by_0 = np.divide(orig_spec_0["abs_flux"],realzn_spec_0_0["abs_flux"])
    div_spec_0_by_1 = np.divide(orig_spec_0["abs_flux"],realzn_spec_0_1["abs_flux"])
    #div_spec_1_by_0 = np.divide(orig_spec_1["abs_flux"],realzn_spec_1_0["abs_flux"])
    #div_spec_1_by_1 = np.divide(orig_spec_1["abs_flux"],realzn_spec_1_1["abs_flux"])

    # are the original and realization spectra really of the same amplitude?
    assert round(np.median(div_spec_0_by_0), 1) == 1.0
    assert round(np.median(div_spec_0_by_1), 1) == 1.0
    #assert round(np.median(div_spec_1_by_0), 1) == 1.0
    #assert round(np.median(div_spec_1_by_1), 1) == 1.0

    ## ## noise injection level not well tested yet

def test_create_norm_spec():
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

    # create noise-churned realizations for each spectrum
    test_input_name_list = [
                            "test_bkgrnd_output_575020m05_noise_ver_000.smo",
                            "test_bkgrnd_output_575020m05_noise_ver_001.smo",
                            ]


    #test_input_name_list = config_red["data_dirs"]["TEST_DIR_SRC"] + "test_input_file_list.list"
    ## ## note TEST_BIN might be changed to TEST_DIR_REZNS_SPEC_NORM, if it can be changed across multiple functions
    test_new_name_list = create_spec_realizations.create_norm_spec(
                            name_list=test_input_name_list,
                            normdir=config_red["data_dirs"]["TEST_DIR_REZNS_SPEC"],
                            finaldir=config_red["data_dirs"]["TEST_DIR_REZNS_SPEC_NORM"])

    # new file name list should have full path and same basename
    assert [i.split("/")[-1] for i in test_new_name_list] == test_input_name_list


def test_read_list():

    abs_stem_src = config_red["data_dirs"]["TEST_DIR_SRC"]

    file_name_test = abs_stem_src + "test_input_file_list.list"

    returned_list = create_spec_realizations.read_list(input_list=file_name_test)

    # is a numpy array with the expected first element returned?
    assert isinstance(returned_list, np.ndarray)
    assert returned_list[0] == "575020m05.smo"


def test_write_bckgrnd_input():
    '''
    Test for function to create input file for the bckgrnd program

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

    indir_test = config_red["data_dirs"]["TEST_DIR_SRC"]
    normdir_test = config_red["data_dirs"]["TEST_DIR_BIN"]
    name_list_test = ["dummy_spec_001.dat","dummy_spec_002.dat","dummy_spec_003.dat"]

    bgrnd_input_filename_test = create_spec_realizations.write_bckgrnd_input(
                                        name_list=name_list_test,
                                        indir=indir_test,
                                        normdir=normdir_test)

    # is the pathname of the written file what we would expect?
    # note this doesn't test whether the file itself is written, but that
    # would be easy to isolate
    assert bgrnd_input_filename_test == indir_test + "bckgrnd_input.txt"


def test_read_bkgrnd_spec():
    # can spectra written by background routine be written in as astropy tables?

    abs_stem_src = config_red["data_dirs"]["TEST_DIR_SRC"]

    # loop over a few files
    # (note they should have 3 columns: 1.) wavelength, 2.) flux, 3.) background flux)
    file_name_array = ["575020m05_noise_ver_000.smo","575020m05_noise_ver_001.smo","575020m10_noise_ver_000.smo"]
    for i in range(0,len(file_name_array)):
        # this is a file that looks like what it should after bkgrnd has done its thing
        # (it should have 3 columns)
        file_name_test = abs_stem_src + file_name_array[i]

        # choose a random spectrum from these four
        returned_table = create_spec_realizations.read_bkgrnd_spec(file_name_test)

        # is returned table a non-empty numpy table?
        assert isinstance(returned_table, astropy.table.table.Table)
        assert len(returned_table) > 0
