'''
Consolidates all the data to make it ready for the MCMC
'''

import glob
import sys
import pickle
import numpy as np
import pandas as pd
from pylab import *
from IPython.display import clear_output
from rrlyrae_metallicity.modules2 import *


def graft_feh(pickle_source_dir=config["data_dirs"]["DIR_PICKLE"],
              stars_and_offsets_info_file=config["file_names"]["RRAB_RRAB_OFFSETS"],
              hk_source_dir=config["data_dirs"]["DIR_SRC"],
              synthetic=False):
    '''
    Read in the EW and phase data, and attach Fe/H values

    INPUTS:
    pickle_source_dir: source directory of pickled mapped Fe/H data
    stars_and_offsets_info_file: file name of star subtypes and their offsets
    hk_source_dir: source directory for HK info
    synthetic:
        =False: the Fe/H will be derived from the basis
        =True: the Fe/H will be extracted from the file name

    OUTPUTS:
    (written pickle file of EW and Fe/H values for each star)
    '''

    ## ## TACK PHASES ONTO LIST OF EWS FROM SPECTRA
    ## ## NEED TO GET RID OF THE 'FAKE' AT SOME POINT

    # read in star names first
    ## ## N.b. this is just the RRabs with RRab offsets for now
    real_data_1 = pickle.load( open(pickle_source_dir  +
                                    stars_and_offsets_info_file, "rb" ) )

    # arrange the data in a way we can use
    # N.b. This is NOT fake data; I'm just appropriating the old variable name
    ## ## Note the ersatz Layden errors for now; need to revisit this with values from his paper
    data_1 = { "star_name": real_data_1[0]["name_star"],
                "feh_lit": real_data_1[0]["feh_highres"],
                "feh_layden": real_data_1[0]["feh_basis"],
                "err_feh_lit": np.zeros(len(real_data_1[0]["feh_basis"])),
                "err_feh_layden": 0.07*np.ones(len(real_data_1[0]["feh_basis"]))}
    #dataset_1 = pd.DataFrame(data=data_1)

    # loop over each star to read in the calculated metallicities
    final_star_feh = pd.DataFrame(columns=["star_name",
                                           "final_feh_center",
                                           "final_feh_lower",
                                           "final_feh_upper"])
    for t in range(0, len(data_1["star_name"])):
        this_star = data_1["star_name"][t]

        # replace space with underscore
        name_star_underscore = str(this_star).replace(" ", "_")

        # read the mapped Fe/H values
        pickle_read_name = (pickle_source_dir + "plot_info_" +
                            name_star_underscore + ".pkl")
        with open(pickle_read_name, 'rb') as f:
                name_star,feh_mapped_array,x_vals,y_vals,xvals_interp,cdf_gauss_info,\
                  idx,idx_1sig_low,idx_1sig_high,shortest_xrange_lower,\
                  shortest_xrange_upper,shortest_xrange_halfway = pickle.load(f)
        this_feh_center = shortest_xrange_halfway
        this_feh_lower = shortest_xrange_lower
        this_feh_upper = shortest_xrange_upper

        final_star_feh = final_star_feh.append({"star_name_underscore": name_star_underscore,
                                               "final_feh_center": this_feh_center,
                                               "final_feh_lower": this_feh_lower,
                                               "final_feh_upper": this_feh_upper},
                                               ignore_index=True)

    # read in the EW and phase info
    hk_ews = pd.read_csv(hk_source_dir + config["file_names"]["MORE_REALISTIC"])

    # paste the feh values onto the HK table
    # loop over each row of the HK table and assign an FeH based
    # on string in empirical spectrum name
    hk_ews["final_feh_center"] = np.nan
    hk_ews["final_feh_lower"] = np.nan
    hk_ews["final_feh_upper"] = np.nan
    hk_ews["Teff"] = np.nan
    hk_ews["logg"] = np.nan

    # loop over each star name (of which our program stars are a subset)
    # and paste the FeH values to the HK table rows corresponding to the
    # empirical spectra for that star
    if not synthetic:
        for star_num in range(0, len(final_star_feh["star_name_underscore"])):
            this_star = final_star_feh["star_name_underscore"][star_num]
            print("Retrieving calculated Fe/H value for " + this_star)
            feh_center_this_star = final_star_feh["final_feh_center"][star_num]
            feh_lower_this_star = final_star_feh["final_feh_lower"][star_num]
            feh_upper_this_star = final_star_feh["final_feh_upper"][star_num]

            # loop over each of our program stars; i.e., empirical spectra
            for em_spec_num in range(0, len(hk_ews["empir_spec_name"])):
                
                # if the star assigned to an FeH value appears in the empirical spectrum name
                if (this_star in hk_ews["empir_spec_name"][em_spec_num]):
                    print("this_star is in hk_ews")
                    hk_ews["final_feh_center"].iloc[em_spec_num] = feh_center_this_star
                    hk_ews["final_feh_lower"].iloc[em_spec_num] = feh_lower_this_star
                    hk_ews["final_feh_upper"].iloc[em_spec_num] = feh_upper_this_star

                #else:
                #    print("Spectrum is not synthetic and does not appear among the program stars.")

    # if these are just synthetic spectra
    elif synthetic:
        print("SYNTHETIC")
        for synth_spec_num in range(0, len(hk_ews["empir_spec_name"])):
            print("Num " + str(synth_spec_num) + " out of " + str(len(hk_ews["empir_spec_name"])))
            this_synth_spectrum_name = hk_ews["empir_spec_name"][synth_spec_num]

            # Fe/H is negative (given RW's naming convention)
            if ("m" in this_synth_spectrum_name):
                feh_center_this_star = -0.1*np.float(this_synth_spectrum_name.split("m")[1])
                feh_lower_this_star = np.float(feh_center_this_star)
                feh_upper_this_star = np.float(feh_center_this_star)
            # Fe/H is positive
            elif ("p" in this_synth_spectrum_name):
                feh_center_this_star = 0.1*np.float(this_synth_spectrum_name.split("p")[1])
                feh_lower_this_star = np.float(feh_center_this_star)
                feh_upper_this_star = np.float(feh_center_this_star)
            Teff_this_star = int(this_synth_spectrum_name.split("m")[0][0:4])
            logg_this_star = np.divide(int(this_synth_spectrum_name.split("m")[0][4:6]),10)

            hk_ews["final_feh_center"].iloc[synth_spec_num] = feh_center_this_star
            hk_ews["final_feh_lower"].iloc[synth_spec_num] = feh_lower_this_star
            hk_ews["final_feh_upper"].iloc[synth_spec_num] = feh_upper_this_star
            hk_ews["Teff"].iloc[synth_spec_num] = Teff_this_star
            hk_ews["logg"].iloc[synth_spec_num] = logg_this_star

        ### START HERE
        # now remove the synthetic spectra with Teff<6000, Teff>7500, logg=2
        # ... this functionality should be moved
        hk_ews = hk_ews[hk_ews.Teff >= 6000]
        hk_ews = hk_ews[hk_ews.Teff <= 7500]
        hk_ews = hk_ews[hk_ews.logg > 2]
        hk_ews = hk_ews.reset_index()
        # ... then find Fe/H abcd; this is 'infinte S/N'
        # ... then check the spectra with xmgrace to see why there are some error large KH bars 
        # ... then put in more realistic noise and repeat

    # fyi
    hk_ews.to_csv('junk.csv')

    print("HK_EWS")
    print(hk_ews)

    # pickle the table of H,K,phases,Fe/H
    ## ## NEED TO ADD STAR TYPE, TOO
    pickle_write_name = pickle_source_dir + config["file_names"]["KH_FINAL_PKL"]
    with open(pickle_write_name, "wb") as f:
        pickle.dump(hk_ews, f)

    return


def winnow_by_phase_and_type(pickle_source_dir=config["data_dirs"]["DIR_PICKLE"],
                         hk_winnowed_write_dir=config["data_dirs"]["DIR_BIN"],
                         remove_rrl_subtype="c"):
    '''
    This removes the program star spectra which are in the bad phase region,
    or are the wrong RRL subtype

    INPUTS:
    pickle_source_dir: directory containing
    hk_winnowed_write_dir: directory to which csv info on H, K is written
    remove_rrl_subtype: RR Lyrae subtype to remove from analysis ("ab", "c", or none)
    '''

    # read in phase boundaries
    min_good, max_good = phase_regions()

    # restore pickle file with all the H,K data
    hk_data = pickle.load( open( pickle_source_dir +
                                 config["file_names"]["KH_FINAL_PKL"], "rb" ) )
    #hk_data_df = pd.DataFrame(hk_data)
    print(hk_data)
    print("hk_data keys:")
    print(hk_data.keys())
    print(hk_data["phase"])

    # drop bad phases
    ## ## NOTE THAT THE DROPNA HERE SEEMS TO BE DROPPING ALL ROWS WITH ANY NANS IN IT (SOME OF THE RRC FEHS ARE NANS)
    ## ## ALSO CHECK THAT WERE NOT LOSING V535 OR V445 THROUGH SILLY NAME DIFFERENCES
    hk_data_winnowed_phase = hk_data.where(np.logical_and(hk_data["phase"] > min_good,
                                                    hk_data["phase"] < max_good)).dropna().reset_index()

    ## ## drop by type, too?
    #hk_data_winnowed_phase_and_type = hk_data_winnowed_phase.where(np.logical_and(hk_data["phase"] > min_good,
    #                                            hk_data["phase"] < max_good)).dropna().reset_index()
    hk_data_winnowed_phase_and_type = hk_data_winnowed_phase


    #hk_data_winnowed_file_name = "hk_data_winnowed.csv"
    hk_data_winnowed_phase.to_csv(hk_winnowed_write_dir +
                                  config["file_names"]["KH_WINNOWED_PHASE_ONLY_FILE_NAME"])
    hk_data_winnowed_phase_and_type.to_csv(hk_winnowed_write_dir +
                                  config["file_names"]["KH_WINNOWED_PHASE_SUBTYPE_FILE_NAME"])

    

    ## ## NEED TO WINNOW BY STAR TYPE, TOO

    return
