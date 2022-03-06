'''
Reads in literature metallicities and makes new Fe/H basis
'''

import pickle
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astroquery.simbad import Simbad
from . import *

class LitFehRaw():
    '''
    Read in Fe/H values from the literature, before making any transformations
    '''

    def __init__(self):
        # map the raw data to object

        # source_dir=config_red["data_dirs"]["DIR_LIT_HIGH_RES_FEH"]):
        source_dir = "/Users/bandari/Documents/git.repos/rrlfe/src/high_res_feh/"

        # stand-in that consists of our program star names
        self.df_our_program_stars = pd.read_csv(source_dir + "our_program_stars_names_only.csv")

        # Fe/H from Layden+ 1994; this may serve as the common basis for RRabs
        self.df_layden_feh = pd.read_csv(source_dir + "layden_1994_abundances.dat")
        # RES: "rather low"

        # Fe/H Clementini+ 1995
        self.df_clementini_feh = pd.read_csv(source_dir + "clementini_1995_abundances.dat")

        # Fe/H Fernley+ 1996
        self.df_fernley96_feh = pd.read_csv(source_dir + "fernley_1996_abundances.dat")
        # RES: 60,000, FeI & FeII, 5900-8100 A

        # Fe/H from Fernley+ 1997
        self.df_fernley97_feh = pd.read_csv(source_dir + "fernley_1997_abundances.dat")
        # RES: 60,000, two FeII lines, 5900-8100 A

        # log(eps) from Lambert+ 1996
        self.df_lambert_logeps = pd.read_csv(source_dir + "lambert_1996_abundances.dat")
        # RES: ~23,000, FeII + photometric models, 3600-9000 A

        # Fe/H from Wallerstein and Huang 2010, arXiv 1004.2017
        self.df_wallerstein_feh = pd.read_csv(source_dir + "wallerstein_huang_2010_abundances.dat")
        # RES: ~30,000, FeII

        # Fe/H from Chadid+ 2017 ApJ 835.2:187 (FeI and II lines)
        self.df_chadid_feh = pd.read_csv(source_dir + "chadid_2017_abundances.dat")
        # RES: 38000, FeI & FeII, 3400-9900 A

        # Fe/H from Liu+ 2013 Res Ast Astroph 13:1307
        self.df_liu_feh = pd.read_csv(source_dir + "liu_2013_abundances.dat")
        # RES: ~60,000, FeI (& FeII?), 5100-6400 A

        # Fe/H from Nemec+ 2013
        self.df_nemec_feh = pd.read_csv(source_dir + "nemec_2013_abundances.dat")
        # RES: ~65,000 or 36,000, FeI & FeII, 5150-5200 A

        # Fe/H from Solano+ 1997
        self.df_solano_feh = pd.read_csv(source_dir + "solano_1997_abundances.dat")
        # RES: 22,000 & 19,000, strong FeI lines, 4160-4390 & 4070-4490 A

        # Fe/H from Pancino+ 2015 MNRAS 447:2404
        self.df_pancino_feh = pd.read_csv(source_dir + "pancino_2015_abundances.dat")
        # RES: >30,000, FeI (weighted average), 4000-8500 A

        # Fe/H from Sneden+ 2017
        self.df_sneden_feh = pd.read_csv(source_dir + "sneden_2017_abundances.dat", delimiter="|")
        # RES: ~27,000 (at 5000 A), FeI & FeII, 3400-9000 A

        # Fe/H from Kemper+ 1982; this might serve as the common basis for RRcs
        self.df_kemper_feh = pd.read_csv(source_dir + "kemper_1982_abundances.dat")

        # Fe/H from Govea+ 2014
        ## ## note: Govea+ has abundances for each phase value, and this
        ## ## includes NLTE phases; how to get single Fe/H?
        self.df_govea_feh = pd.read_csv(source_dir + "govea_2014_abundances.dat")

def map_names(df_pass):
    # find common ASAS names

    # convert output in astropy.Table format to DataFrame
    df_test = pd.DataFrame(test.df_govea_feh)

    # make column of long-form, ASAS designations
    df_pass["ASAS_common"] = np.nan

    # loop over rows, parse as necessary
    for row_num in range(0,len(df_pass)):
        name_initial = df_pass["name"]




def main():

    # read in raw
    test = LitFehRaw()
    print(test)
    import ipdb; ipdb.set_trace()

    # make transformations to get single Fe/H value

    # expand abbreviated ASAS names (ex. Govea)

    # get common ASAS names from Simbad, put them in new col
    # (note that will have to remove '_' in many strings)
    '''
    .values_equal()

    # make common basis
    sdiff = astropy.table.join(table_1, table_2, keys=["ASAS_common"], join_type="inner")
    '''


# entry point
if __name__ == '__main__':
    sys.exit(main())
