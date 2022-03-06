import matplotlib
matplotlib.use('Agg')

import sys, os
from configparser import ConfigParser, ExtendedInterpolation
import pandas as pd
import astropy

current_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(current_dir, "../"))
sys.path.insert(0, target_dir)

# import more things with changed system path
from modules import *
from modules import scrape_ew_and_errew
from conf import *
import numpy as np
import glob

# configuration data for reduction
config_red = ConfigParser(interpolation=ExtendedInterpolation()) # for parsing values in .init file
# config for reduction to find a, b, c, d
config_red.read(os.path.join(os.path.dirname(__file__), '../conf', 'config_red.ini'))


def test_line_order_check():

    dummy_lines = [3933.660,3970.075,4101.7100,4340.472,4861.290]
    test_num_glitches_0 = scrape_ew_and_errew.line_order_check(dummy_lines)
    dummy_lines[1] = 3990.
    dummy_lines[3] = 4320.
    test_num_glitches_1 = scrape_ew_and_errew.line_order_check(dummy_lines)
    dummy_lines[0] = 3913.
    dummy_lines[2] = 4121.
    dummy_lines[4] = 4881.
    test_num_glitches_2 = scrape_ew_and_errew.line_order_check(dummy_lines)

    # assert glitches is boolean
    assert test_num_glitches_0 == 0
    assert test_num_glitches_1 == 1
    assert test_num_glitches_2 == 1


def test_Scraper():

    '''
    write_dir_test = config_red["data_dirs"]["TEST_DIR_BIN"]
    robo_dir = config_red["sys_dirs"]["DIR_ROBO"]
    file_names_test = glob.glob(config_red["data_dirs"]["TEST_DIR_SRC"] + "spec_norm_final/*")
    '''

    # instantiate
    scraper_instance = scrape_ew_and_errew.Scraper(subdir = config_red["data_dirs"]["TEST_DIR_ROBO_OUTPUT"],
                                                   file_scraped_info = config_red["data_dirs"]["TEST_DIR_BIN"]+"scraper_output/"+config_red["file_names"]["SCRAPED_EW_ALL_DATA"])


    # try a single instance; does it work?
    # note the writing of files is not directly tested here
    function_state = True
    try:
        test = scraper_instance()
    except Exception as e:
        # e contains printable attributes of exception object
        function_state = False

    # assert that instantiation worked
    assert function_state

    # make sure lines are really being identified correctly
    assert np.allclose(test.where(test["line_name"]=="CaIIK").dropna()["wavel_found_center"],3933.660, atol=2.)
    assert np.allclose(test.where(test["line_name"]=="Heps").dropna()["wavel_found_center"],3970.075, atol=2.)
    assert np.allclose(test.where(test["line_name"]=="Hdel").dropna()["wavel_found_center"],4101.710, atol=2.)
    assert np.allclose(test.where(test["line_name"]=="Hgam").dropna()["wavel_found_center"],4340.472, atol=2.)
    assert np.allclose(test.where(test["line_name"]=="Hbet").dropna()["wavel_found_center"],4861.290, atol=2.)

    # only 2 of the 3 spectra should have been scraped, because one should have
    # triggered a parsing errors
    assert len(test.where(test["line_name"]=="CaIIK").dropna()) == 2



def test_quality_check():

    data_out_test = scrape_ew_and_errew.quality_check(
                        read_in_filename = config_red["data_dirs"]["TEST_DIR_BIN"]+"scraper_output/"+config_red["file_names"]["TEST_SCRAPED_EW_ALL_DATA"],
                        write_out_filename = config_red["data_dirs"]["TEST_DIR_BIN"]+"scraper_output/"+config_red["file_names"]["TEST_SCRAPED_EW_DATA_GOOD_ONLY"])

    # lots of checks of data types
    # note this uses .iloc[0] instead of [0], because bad rows with index 0 may
    # have been removed
    assert isinstance(data_out_test["wavel_stated_center"].iloc[0],np.float64)
    assert isinstance(data_out_test["wavel_found_center"].iloc[0],np.float64)
    assert isinstance(data_out_test["gaussianSigma"].iloc[0],np.float64)
    assert isinstance(data_out_test["gaussianAmp"].iloc[0],np.float64)
    assert isinstance(data_out_test["uncertaintyMu"].iloc[0],np.float64)
    assert isinstance(data_out_test["uncertaintySigma"].iloc[0],np.float64)
    assert isinstance(data_out_test["uncertaintyAmp"].iloc[0],np.float64)
    assert isinstance(data_out_test["priorMu"].iloc[0],np.float64)
    assert isinstance(data_out_test["priorSigma"].iloc[0],np.float64)
    assert isinstance(data_out_test["priorAmp"].iloc[0],np.float64)
    assert isinstance(data_out_test["EQW"].iloc[0],np.float64)
    assert isinstance(data_out_test["uncertaintyEQW"].iloc[0],np.float64)
    assert isinstance(data_out_test["chiSqr"].iloc[0],np.float64)
    assert isinstance(data_out_test["flags"].iloc[0],str)
    assert isinstance(data_out_test["blendGroup"].iloc[0],np.int64)
    assert isinstance(data_out_test["line_name"].iloc[0],str)
    assert isinstance(data_out_test["robolines_file_name"].iloc[0],str)
    assert isinstance(data_out_test["realization_spec_file_name"].iloc[0],str)
    assert isinstance(data_out_test["quality"].iloc[0],str)

def test_stack_spectra():

    print("input list")
    print(config_red["data_dirs"]["TEST_DIR_SRC"] + "test_input_file_list.list")
    print("read in file name")

    data_stacked_test = scrape_ew_and_errew.stack_spectra(
                            read_in_filename = config_red["data_dirs"]["TEST_DIR_BIN"]+"scraper_output/"+config_red["file_names"]["TEST_SCRAPED_EW_DATA_GOOD_ONLY"],
                            write_out_filename = config_red["data_dirs"]["TEST_DIR_BIN"]+"scraper_output/"+config_red["file_names"]["TEST_RESTACKED_EW_DATA_GOOD_ONLY"],
                            input_list = config_red["data_dirs"]["TEST_DIR_SRC"] + "test_input_file_list.list")

    print("data_stacked")
    print(data_stacked_test.keys())

    # lots of checks of data types
    # note this uses .iloc[0] instead of [0], because bad rows with index 0 may
    # have been removed

    assert isinstance(data_stacked_test["realization_spec_file_name"].iloc[0],str)
    assert isinstance(data_stacked_test["orig_spec_file_name"].iloc[0],str)
    assert isinstance(data_stacked_test["EW_Hbeta"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["err_EW_Hbeta_from_robo"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["EW_Hdelta"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["err_EW_Hdelta_from_robo"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["EW_Hgamma"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["err_EW_Hgamma_from_robo"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["EW_Heps"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["err_EW_Heps_from_robo"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["EW_CaIIK"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["err_EW_CaIIK_from_robo"].iloc[0],np.float64)


def test_generate_net_balmer():

    # generate the fake data: H_del =

    params_data, data_net_balmer_test = scrape_ew_and_errew.generate_net_balmer(
                                                    read_in_filename = config_red["data_dirs"]["TEST_DIR_BIN"]+"scraper_output/test_stacked_data_pre_net_balmer_calc.csv",
                                                    write_out_filename = config_red["data_dirs"]["TEST_DIR_BIN"]+"scraper_output/test_stacked_data_post_net_balmer_calc.csv")

    # is the Balmer line a true element wise average?
    print("-----")
    print(data_net_balmer_test["EW_Hgamma"])
    print(data_net_balmer_test["EW_Hdelta"])
    print(data_net_balmer_test["EW_Balmer"])

    print(type(np.mean([data_net_balmer_test["EW_Hgamma"],data_net_balmer_test["EW_Hdelta"]], axis=0)))
    print(np.array(data_net_balmer_test["EW_Balmer"]))
    print(np.testing.assert_allclose(np.array(data_net_balmer_test["EW_Balmer"]),np.mean([data_net_balmer_test["EW_Hgamma"],data_net_balmer_test["EW_Hdelta"]], axis=0)))

    #assert np.array(data_net_balmer_test["EW_Balmer"]) == np.mean([data_net_balmer_test["EW_Hgamma"],data_net_balmer_test["EW_Hdelta"]], axis=0)
    '''
    # check data type of newly-added data
    assert isinstance(data_net_balmer_test["EW_Balmer"].iloc[0],np.float64)
    assert isinstance(data_net_balmer_test["err_EW_Balmer_based_Robo"].iloc[0],np.float64)
    '''


def test_generate_addl_ew_errors():
    # placeholder for now, until more decisions about how to calculate EW errors

    test_df_postbalmer_errors = scrape_ew_and_errew.generate_addl_ew_errors(read_in_filename=config_red["data_dirs"]["TEST_DIR_SRC"]+config_red["file_names"]["TEST_RESTACKED_EW_DATA_W_NET_BALMER"],
                                                        write_out_filename=config_red["data_dirs"]["TEST_DIR_BIN"]+config_red["file_names"]["TEST_RESTACKED_EW_DATA_W_NET_BALMER_ERRORS"])

    # loop through batches of rows corresponding to an individual spectrum, and
    # make sure the errors are consistent and the value expected
    array_1 = test_df_postbalmer_errors["err_EW_Balmer_based_noise_churning"].where(test_df_postbalmer_errors["orig_spec_file_name"]=="575020m10.smo").dropna().values
    array_2 = test_df_postbalmer_errors["err_EW_Balmer_based_noise_churning"].where(test_df_postbalmer_errors["orig_spec_file_name"]=="575020m15.smo").dropna().values
    array_3 = test_df_postbalmer_errors["err_EW_Balmer_based_noise_churning"].where(test_df_postbalmer_errors["orig_spec_file_name"]=="575020m20.smo").dropna().values
    # check that error values are same
    assert np.all(array_1)
    assert np.all(array_2)
    assert np.all(array_3)
    # check that one value has expected value
    assert round(array_1[0], 3) == 0.023
    assert round(array_2[0], 3) == 0.020
    assert round(array_3[0], 3) == 0.048


def test_add_synthetic_meta_data():

    combined_data = scrape_ew_and_errew.add_synthetic_meta_data(input_list = config_red["data_dirs"]["TEST_DIR_SRC"] + config_red["file_names"]["TEST_LIST_SPEC_PHASE"],
                                                read_in_filename = config_red["data_dirs"]["TEST_DIR_EW_PRODS"]+config_red["file_names"]["TEST_RESTACKED_EW_DATA_W_NET_BALMER_ERRORS"],
                                                write_out_filename = config_red["data_dirs"]["TEST_DIR_EW_PRODS"]+config_red["file_names"]["TEST_RESTACKED_EW_DATA_W_METADATA_WRITEOUT"])

    # columns from Robospect output and meta-data are all there
    assert "wavel_stated_center" in combined_data.columns
    assert "feh" in combined_data.columns
    assert "teff" in combined_data.columns
    # there are no nans in the table
    assert np.sum(combined_data.isnull().sum()) == 0
