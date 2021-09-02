import matplotlib
matplotlib.use('Agg')

import sys, os
import configparser
import pandas as pd
import astropy

current_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(current_dir, "../"))
sys.path.insert(0, target_dir)

# import more things with changed system path
from modules import *
from modules import run_robo
from conf import *
import numpy as np
import glob

# configuration data for reduction
config_red = configparser.ConfigParser() # for parsing values in .init file
# config for reduction to find a, b, c, d
config_red.read(os.path.join(os.path.dirname(__file__), '../conf', 'config_red.ini'))


def test_RunRobo():

    write_dir_test = config_red["data_dirs"]["TEST_DIR_BIN"]
    robo_dir = config_red["data_dirs"]["DIR_ROBO"]
    file_names_test = glob.glob(config_red["data_dirs"]["TEST_DIR_SRC"] + "spec_norm_final/*")

    # instantiate
    run_robospect_instance = run_robo.RunRobo(write_dir = write_dir_test, robo_dir = robo_dir)

    # try a single instance; does it work?
    # note the writing of files is not directly tested here
    function_state = True
    try:
        run_robospect_instance(file_names_test[0])
    except Exception as e:
        # e contains printable attributes of exception object
        function_state = False

    assert function_state