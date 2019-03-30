import configparser
import os

# configuration data
config = configparser.ConfigParser() # for parsing values in .init file
config.read("rrlyrae_metallicity/modules2/config.ini")

from setuptools import Distribution
from setuptools.command.install import install

# The class OnlyGetScriptPath() and function get_setuptools_script_dir() are from the setup.py script
# in the Apogee repository by jobovy
# https://github.com/jobovy/apogee/blob/master/setup.py

class OnlyGetScriptPath(install):
    def run(self):
        self.distribution.install_scripts = self.install_scripts

def get_setuptools_script_dir():
    " Get the directory setuptools installs scripts to for current python "
    dist = Distribution({'cmdclass': {'install': OnlyGetScriptPath}})
    dist.dry_run = True  # not sure if necessary
    dist.parse_config_files()
    command = dist.get_command_obj('install')
    command.ensure_finalized()
    command.run()
    return dist.install_scripts

def make_dirs():
    '''
    Make directories for housing files/info if they don't already exist
    '''

    # loop over all directory paths we will need
    for vals in config["data_dirs"]:
        abs_path_name = str(config["data_dirs"][vals])
        print("Directory exists: " + abs_path_name)
        
        # if directory does not exist, create it
        if not os.path.exists(abs_path_name):
            os.makedirs(abs_path_name)
            print("Made directory " + abs_path_name)

def phase_regions():
    '''
    Read in the boundary between good and bad phase regions
    '''

    # obtain values as floats
    value1 = config.getfloat("phase", "MIN_GOOD")
    value2 = config.getfloat("phase", "MAX_GOOD")
    
    return value1, value2
