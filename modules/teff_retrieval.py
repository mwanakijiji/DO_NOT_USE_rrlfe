#!/usr/bin/env python
# coding: utf-8

# This makes plots showing the effective temperature retrievals based on synthetic spectra
# produced by R.W.

# Created from parent restacking_scraped_data.ipynb 2021 March 17 by E.S.

import pandas as pd
from astropy.io.fits import getdata
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os.path
import csv
import git
from . import *


def line_fit_temp_range(x_data_pass, y_data_pass, t_min, t_max):
    '''
    Find line of best fit

    INPUTS:
    x_data_pass: abcissa
    y_data_pass: ordinate

    OUTPUTS:
    m:      slope
    err_m:  error in slope
    b:      y-intercept
    err_b:  error in y-intercept
    '''

    # remove the stuff outside of 6000-7250 K
    #x_data_rrl = x_data_pass.where(np.logical_and(x_data_pass>=5900,x_data_pass<=7350))
    #y_data_rrl = x_data_pass.where(np.logical_and(x_data_pass>=5900,x_data_pass<=7350))
    x_data_rrl = x_data_pass[np.where(np.logical_and(y_data_pass>=t_min,y_data_pass<=t_max))]
    y_data_rrl = y_data_pass[np.where(np.logical_and(y_data_pass>=t_min,y_data_pass<=t_max))]

    coeff, cov = np.polyfit(x_data_rrl, y_data_rrl, 1, full=False, cov=True)
    m = coeff[0]
    b = coeff[1]
    err_m = np.sqrt(np.diag(cov))[0]
    err_b = np.sqrt(np.diag(cov))[1]

    logging.info("Fitting a Teff vs. Balmer line trend. Temperature range "+\
                    "restricted to " + str(int(t_min)) + ", " + str(int(t_max)) +" K")

    return m, err_m, b, err_b


def temp_vs_balmer(df_poststack_file_name_read = config_red["data_dirs"]["DIR_EW_PRODS"]+config_red["file_names"]["RESTACKED_EW_DATA_W_METADATA"],
                    df_poststack_file_name_write = config_red["data_dirs"]["DIR_EW_PRODS"] + config_red["file_names"]["RESTACKED_EW_DATA_GOOD_ONLY_TEFFFIT"],
                    plot_write = config_red["data_dirs"]["DIR_BIN"] + config_red["file_names"]["PLOT_TEFF_VS_BALMER"],
                    teff_data_write = config_red["data_dirs"]["DIR_BIN"] + config_red["file_names"]["TREND_TEFF_VS_BALMER"],
                    plot = True,
                    test_flag=False):
    '''
    Finds a linear Teff vs. Balmer EW relation. This is an ancillary step before
    running the MCMC further downstream in the pipeline.

    INPUTS:
    df_poststack_file_name_read: name of file that contains all the data from the upstream
        pipeline and will be read in for the fit; it should contain columns with 'teff'
        and 'EW_Balmer', with which a simple linear fit is made
    df_poststack_file_name_write: name of file to write; this file is the same as
        the one read in, except that now it also includes the best-fit values of the Teff
    teff_data_write: file name of txt file containing info on the lienar trend
    plot_write: file name of Teff vs Balmer plot to write
    teff_data_write: name of file to write Teff data to
    plot: flag whether to write plot or not
    test_flag: if testing, suppress prompts in terminal

    OUTPUTS:
    m:      slope
    err_m:  error in slope
    b:      y-intercept
    err_b:  error in y-intercept
    '''

    # the min and max Teff of spectra that the linear fit will be made to
    t_min = int(config_red["teff_linear"]["MIN_TEFF"])
    t_max = int(config_red["teff_linear"]["MAX_TEFF"])

    # read in data
    df_poststack = pd.read_csv(df_poststack_file_name_read)

    # find linear trend of net Balmer EW with Teff
    teff = df_poststack["teff"].values.astype(float)
    # fit a straight line: net Balmer
    ews_Balmer = df_poststack["EW_Balmer"].values.astype(float)

    m, err_m, b, err_b = line_fit_temp_range(x_data_pass=ews_Balmer,
                                                y_data_pass=teff,
                                                t_min=t_min,
                                                t_max=t_max)

    logging.info("Best-fit line for Teff=m*EW_Balmer + b is [m, err_m, b, err_b] = " +
                "[" + str(m) + ", " + str(err_m) + ", " + str(b) + ", " + str(err_b) + "]")

    # add the best-fit Teffs to dataframe
    teffs_bestfit = np.add(np.multiply(m,ews_Balmer),b)
    df_poststack["teff_bestfit"] = teffs_bestfit

    # write all the data to file
    df_poststack.to_csv(df_poststack_file_name_write,index=False)
    logging.info("Wrote out data file including linear-best-fit Teffs to " + df_poststack_file_name_write)

    # retrieve hash
    repo = git.Repo(search_parent_directories=True)
    sha = repo.head.object.hexsha

    # arrange info into dictionary

    linfit_info={
                "Hash"      :   sha,
                "Teff_min"  :   t_min,
                "Teff_max"  :   t_max,
                "m"         :      m,
                "err_m"     :  err_m,
                "b"         :      b,
                "err_b"     :  err_b
                }

    # write the Teff trend parameters alone to a separate text file
    # (note this overwrites any previous existing file)
    if (os.path.exists(teff_data_write) and test_flag==False): # pragma: no cover

        print(teff_data_write)
        input("Text file containing Teff linear fit trend already exists! \n" + \
                teff_data_write + "\n" + \
                "Do what you want with that file and hit [ENTER] (will overwrite)")

    with open(teff_data_write, 'w') as file1:

        # header
        file1.write("Linear fit to Teff vs Balmer EW\n")

        # data
        for key, value in linfit_info.items():
            file1.write('%s:%s\n' % (key, value))

    file1.close()
    logging.info("Wrote out text file with linear-best-fit params to " + teff_data_write)

    # save an FYI plot
    if plot: # pragma: no cover
        plt.clf()
        plt.title("Teff from the Balmer EW\n[m, err_m, b, err_b] = \n" +
                    "[" + str(np.round(m,2)) + ", " + str(np.round(err_m,2)) + ", " + str(np.round(b,2)) + ", " + str(np.round(err_b,2)) + "]")
        plt.axhline(y=t_max, color="k", alpha=0.4)
        plt.axhline(y=t_min, color="k", alpha=0.4)
        plt.plot(ews_Balmer, teffs_bestfit, linestyle='--')
        plt.scatter(ews_Balmer, teff, color="k", s=3)

        plt.ylabel("Teff (K)")
        plt.xlabel("EW (Angstr)")
        plt.tight_layout()
        plt.savefig(plot_write)
        plt.clf()
        logging.info("Wrote out plot of Teffs vs. Balmer EW to " + plot_write)

    return df_poststack
