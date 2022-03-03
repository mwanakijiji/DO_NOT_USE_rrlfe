'''
Initialization
'''

from configparser import ConfigParser, ExtendedInterpolation
import os
import sys
import multiprocessing
import logging
from setuptools import Distribution
from setuptools.command.install import install
from datetime import datetime

# set up logging, to print to screen and save to file simultaneously
timestring_start = datetime.now().strftime("%Y%m%d_%H%M%S")
log_filename = timestring_start + ".log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s.%(msecs)03d [%(levelname)s] %(module)s - %(funcName)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(log_filename),
        logging.StreamHandler(sys.stdout)
    ]
)

# configuration data for reduction
config_red = ConfigParser(interpolation=ExtendedInterpolation()) # for parsing values in .init file
# config for reduction to find a, b, c, d
config_red.read(os.path.join(os.path.dirname(__file__), '../conf', 'config_red.ini'))

# config for applying a calibration
config_apply = ConfigParser(interpolation=ExtendedInterpolation())

config_apply.read(os.path.join(os.path.dirname(__file__), '../conf', 'config_apply.ini'))

# set pathnames for important files that are used by different modules
cc_bkgrnd_file_path_abs = config_red["data_dirs"]["DIR_SRC"] + "/bkgrnd.cc"
compiled_bkgrnd_file_path_abs = config_red["data_dirs"]["DIR_BIN"] + "/bkgrnd"

ncpu = multiprocessing.cpu_count()

# set some constants

# vestigial constant in MCMC (this is NOT an astrophysical Teff)
Teff = 0.0586758
# coefficients from first line of Table 8 in Layden+ 1994
# (reddening not included), to serve as MCMC starting point
a_layden = 13.542
b_layden = -1.072
c_layden = 3.971
d_layden = -0.271


# The class OnlyGetScriptPath() and function get_setuptools_script_dir()
# are from the setup.py script in the Apogee repository by jobovy
# https://github.com/jobovy/apogee/blob/master/setup.py

class OnlyGetScriptPath(install):
    def run(self):
        self.distribution.install_scripts = self.install_scripts

def get_setuptools_script_dir():
    '''
    Get the directory setuptools installs scripts to for current python
    '''
    dist = Distribution({'cmdclass': {'install': OnlyGetScriptPath}})
    dist.dry_run = True  # not sure if necessary
    dist.parse_config_files()
    command = dist.get_command_obj('install')
    command.ensure_finalized()
    command.run()
    return dist.install_scripts

def make_dirs(objective="apply_calib"):
    '''
    Make directories for housing files/info if they don't already exist
    '''
    logging.info("## Making directories ##")

    # make directories for
    # 1. reduction of spectra to find a, b, c, d (objective = "find_calib"), or
    # 2. to apply the solution (objective = "apply_calib"; default)
    print("----------")
    print(objective)
    if (objective == "apply_calib"):
        config_choice = config_apply
    elif (objective == "find_calib"):
        config_choice = config_red

    # loop over all directory paths we will need
    for vals in config_choice["data_dirs"]:
        abs_path_name = str(config_choice["data_dirs"][vals])
        logging.info("Directory exists: " + abs_path_name)

        # if directory does not exist, create it
        if not os.path.exists(abs_path_name):
            original_umask = os.umask(0) # original system permission
            os.makedirs(abs_path_name, 0o777)
            #os.mkdir(abs_path_name)
            logging.info("Made directory " + abs_path_name)
            os.umask(original_umask) # revert to previous permission status

def phase_regions():
    '''
    Read in the boundary between good and bad phase regions
    '''

    # obtain values as floats
    value1 = config_red.getfloat("phase", "MIN_GOOD")
    value2 = config_red.getfloat("phase", "MAX_GOOD")

    return value1, value2
