from rrlyrae_metallicity.modules2 import *
from rrlyrae_metallicity.modules2 import graft_phases
from rrlyrae_metallicity.modules2 import edification

#import sys, os
#myPath = os.path.dirname(os.path.abspath(__file__))
#sys.path.insert(0, myPath + '/../')

#from rrlyrae_metallicity.modules2 import yada_compile_normalization
#from junktest_compile_normalizations import compile_bkgrnd


def test_config():
    assert edification.graft_feh() == '22'

'''
def test_getcwd():
    #junktest_compile_normalizations.compile_bkgrnd()
    # configuration data
    graft_phases.yada_yada() == '22'
    assert True

# check if the directory-making function works
def test_make_dirs():

    # call function to make directories
    make_dirs()

    # do all the directories exist now?
    for vals in config["data_dirs"]:
        abs_path_name = str(config["data_dirs"][vals])
        assert os.path.exists(abs_path_name)

# test if the phase region boundaries are being read in correctly
def test_phase_regions():

    min_good_phase, max_good_phase = phase_regions()

    # is min smaller than max
    assert min_good_phase < max_good_phase

    # are the phases interpreted as floats
    assert isinstance(min_good_phase,float)
'''