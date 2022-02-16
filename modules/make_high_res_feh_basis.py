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

def lit_feh_raw():
    '''
    Read in Fe/H values from the literature, before making any transformations
    '''

    # source_dir=config_red["data_dirs"]["DIR_LIT_HIGH_RES_FEH"]):
    source_dir = "/Users/bandari/Documents/git.repos/rrlfe/src/high_res_feh/"

    # stand-in that consists of our program star names
    df_our_program_stars = pd.read_csv(source_dir + "our_program_stars_names_only.csv")
    print(df_our_program_stars)

    # Fe/H from Layden+ 1994; this may serve as the common basis for RRabs
    df_layden_feh = pd.read_csv(source_dir + "layden_1994_abundances.dat")
    print(df_layden_feh)
    # RES: "rather low"

    # Fe/H Clementini+ 1995
    df_clementini_feh = pd.read_csv(source_dir + "clementini_1995_abundances.dat")
    print(df_clementini_feh)

    # Fe/H Fernley+ 1996
    df_fernley96_feh = pd.read_csv(source_dir + "fernley_1996_abundances.dat")
    print(df_fernley96_feh)
    # RES: 60,000, FeI & FeII, 5900-8100 A

    # Fe/H from Fernley+ 1997
    df_fernley97_feh = pd.read_csv(source_dir + "fernley_1997_abundances.dat")
    print(df_fernley97_feh)
    # RES: 60,000, two FeII lines, 5900-8100 A

    # log(eps) from Lambert+ 1996
    df_lambert_logeps = pd.read_csv(source_dir + "lambert_1996_abundances.dat")
    print(df_lambert_logeps)
    # RES: ~23,000, FeII + photometric models, 3600-9000 A

    # Fe/H from Wallerstein and Huang 2010, arXiv 1004.2017
    df_wallerstein_feh = pd.read_csv(source_dir + "wallerstein_huang_2010_abundances.dat")
    print(df_wallerstein_feh)
    # RES: ~30,000, FeII

    # Fe/H from Chadid+ 2017 ApJ 835.2:187 (FeI and II lines)
    df_chadid_feh = pd.read_csv(source_dir + "chadid_2017_abundances.dat")
    print(df_chadid_feh)
    # RES: 38000, FeI & FeII, 3400-9900 A

    # Fe/H from Liu+ 2013 Res Ast Astroph 13:1307
    df_liu_feh = pd.read_csv(source_dir + "liu_2013_abundances.dat")
    print(df_liu_feh)
    # RES: ~60,000, FeI (& FeII?), 5100-6400 A

    # Fe/H from Nemec+ 2013
    df_nemec_feh = pd.read_csv(source_dir + "nemec_2013_abundances.dat")
    print(df_nemec_feh)
    # RES: ~65,000 or 36,000, FeI & FeII, 5150-5200 A

    # Fe/H from Solano+ 1997
    df_solano_feh = pd.read_csv(source_dir + "solano_1997_abundances.dat")
    print(df_solano_feh)
    # RES: 22,000 & 19,000, strong FeI lines, 4160-4390 & 4070-4490 A

    # Fe/H from Pancino+ 2015 MNRAS 447:2404
    df_pancino_feh = pd.read_csv(source_dir + "pancino_2015_abundances.dat")
    print(df_pancino_feh)
    # RES: >30,000, FeI (weighted average), 4000-8500 A

    # Fe/H from Sneden+ 2017
    df_sneden_feh = pd.read_csv(source_dir + "sneden_2017_abundances.dat", delimiter="|")
    print(df_sneden_feh)
    # RES: ~27,000 (at 5000 A), FeI & FeII, 3400-9000 A

    # Fe/H from Kemper+ 1982; this might serve as the common basis for RRcs
    df_kemper_feh = pd.read_csv(source_dir + "kemper_1982_abundances.dat")
    print(df_kemper_feh)

    # Fe/H from Govea+ 2014
    ## ## note: Govea+ has abundances for each phase value, and this
    ## ## includes NLTE phases; how to get single Fe/H?
    df_govea_feh = pd.read_csv(source_dir + "govea_2014_abundances.dat")
    print(df_govea_feh)

    return df_sneden_feh


def main():

    # read in raw
    test = lit_feh_raw()
    print(test)
    import ipdb; ipdb.set_trace()

    # make transformations to get single Fe/H value

    # expand abbreviated ASAS names (ex. Govea)

    # get common ASAS names from Simbad, put them in new col
    # (note that will have to remove '_' in many strings)

    # make common basis


# entry point
if __name__ == '__main__':
    sys.exit(main())
