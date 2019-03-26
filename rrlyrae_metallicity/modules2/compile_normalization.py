from subprocess import Popen,PIPE,check_call,CalledProcessError
import os
import shutil
from modules2 import *

'''
Compile spectral normalization script
'''

# configuration data
config = configparser.ConfigParser() # for parsing values in .init file
config.read("rrlyrae_metallicity/modules2/config.ini")

def compile_bkgrnd():
    _COMPILE_BKGRND = True
    if _COMPILE_BKGRND:
        if True:
            
            print("Compiling background normalization script...")
            bkgrnd_compile = Popen(["g++","-o",
                                    config["data_dirs"]["DIR_BIN"] + "bkgrnd",
                                    config["data_dirs"]["DIR_SRC"] + "bkgrnd.cc"],
                                    stdout=PIPE,stderr=PIPE)
            
            # now move the executable to the right path (the Python setup tools path)
            # but what if I'm not using a conda environment? I won't know where to move it, beyond just get_setuptools_script_dir()
            shutil.copy(config["data_dirs"]["DIR_BIN"] + "bkgrnd",
                        get_setuptools_script_dir()) 
        try:
            # (why is removal of original file not working without throwing an error??)
            check_call(["rm","-rf",config["data_dirs"]["DIR_BIN"] + "bkgrnd"]) 
        except CalledProcessError:
            print("Binary bkgrnd not removed from original location")
