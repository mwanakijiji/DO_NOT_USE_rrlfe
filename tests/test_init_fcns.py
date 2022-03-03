import matplotlib
matplotlib.use('Agg')

#from ../rrlyrae_metallicity.modules2 import *

import sys, os
#myPath = os.path.dirname(os.path.abspath(__file__))
#sys.path.insert(0, myPath)

current_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(current_dir, "../"))
print(current_dir)
print(target_dir)
sys.path.insert(0, target_dir)

#from rrlyrae_metallicity.rrlyrae_metallicity import modules2
from modules import *
from conf import *
#from rrlyrae_metallicity.rrlyrae_metallicity.modules2 import *


# check if the directory-making function works
def test_make_dirs():

    # call function to make directories
    # (here we choose 'objective' to apply a calibration, though it could
    # also be to find a new calibration)
    make_dirs(objective = "find_calib")

    # do all the directories exist now?
    for vals in config_red["data_dirs"]:
        abs_path_name = str(config_red["data_dirs"][vals])
        assert os.path.exists(abs_path_name)


# test if the phase region boundaries are being read in correctly
def test_phase_regions():

    min_good_phase, max_good_phase = phase_regions()

    # is min smaller than max
    assert min_good_phase < max_good_phase

    # are the phases interpreted as floats
    assert isinstance(min_good_phase,float)
